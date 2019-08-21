/**
 *
 */
package org.theseed.genome.contigs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.sequence.FastaInputStream;
import org.theseed.sequence.Sequence;
import org.theseed.utils.ICommand;

/**
 * This command reads in a FASTA file and produces an input file to request coding information from
 * a model.  The output file is tab-delimited, with the metadata in the first column and the
 * sensors in the remaining columns.
 *
 * The positional parameters are the names of the FASTA files.
 *
 * The following command-line options are supported.
 *
 * -u	the number of positions to examine to the left (upstream) of the target position
 * -d	the number of positions to examine to the right (downstream) of the target position
 * -v	write progress messages to STDERR
 * -f	filter for known edge codons
 *
 * --sensor		type of DNA sensor to use
 * 		direct	each base pair converts to a single number
 * 		codon	each trio of base pairs is converted to a string
 * 		channel	each base pair is converted to a string indicating the base pair
 * 		aminoacid
 * 				each trio of base pairs is converted to its amino acid
 *
 * @author Bruce Parrello
 *
 */
public class FastaProcessor implements ICommand {

    // FIELDS
    /** factory object for creating contig sensors */
    private ContigSensorFactory factory;


    // COMMAND-LINE OPTIONS

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;


    /** sensor width, upstream */
    @Option(name="-u", aliases={"--upstream", "--left"}, metaVar="14", usage="upstream distance for sensors")
    private void setLeftWidth(int newWidth) {
        ContigSensorFactory.setLeftWidth(newWidth);
    }

    /** sensor width, downstream */
    @Option(name="-d", aliases={"--downstream", "--right"}, metaVar="21", usage="downstream distance for sensors")
    private void setRightWidth(int newWidth) {
        ContigSensorFactory.setRightWidth(newWidth);
    }

    /** debug switch */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="write progress messages to STDERR")
    private boolean debug;

    /** filter for edge codons */
    @Option(name="-f", aliases={"--edgeFilter"}, usage="filter for known edge codons")
    private boolean edgeFilter;

    /** sensor type */
    @Option(name="--sensor", metaVar="channel", usage="type of DNA sensor to use")
    private void setFactory(ContigSensorFactory.Type type) {
        this.factory = ContigSensorFactory.create(type);
    }

    /** FASTA file names */
    @Argument(index=0, metaVar="file1.fa file2.fa ...", usage="FASTA files to process")
    private List<File> inFiles;

    /**
     * Parse the command line.
     *
     * @return TRUE if successful, FALSE if the command should be aborted
     */
    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.debug = false;
        this.factory = ContigSensorFactory.create(ContigSensorFactory.Type.CHANNEL);
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // Insure the FASTA files are valid.
                for (File inFile : inFiles) {
                    if (! inFile.exists()) {
                        throw new FileNotFoundException(inFile + " does not exist.");
                    }
                }
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println("Error processing input files: " + e.getMessage());
        }
        return retVal;
    }

    @Override
    public void run() {
        try {
            // Create the output header.  The first column is the metadata location, the second
            // is the codon itself (also metadata) and the remaining columns are sensors.
            System.out.println("Location\tCodon\t" + this.factory.sensor_headers());
            // Set up the codon filter.
            CodonFilter filter = null;
            if (this.edgeFilter)
                filter = new CodonFilter("ATG", "GTG", "TTG", "TAA", "TAG", "TGA");
            // Now we loop through the sequences, producing output.
            for (File inFile : this.inFiles) {
                if (debug) System.err.println("Processing file " + inFile + ".");
                FastaInputStream inStream = new FastaInputStream(inFile);
                for (Sequence inSeq : inStream) {
                    int limit = inSeq.length();
                    String sequence = inSeq.getSequence();
                    // For this sequence, output all the sensors.
                    for (int pos = 1; pos <= limit; pos++) {
                        if (filter == null || filter.matches(pos, sequence)) {
                            ContigSensor sensor = this.factory.create(inSeq.getLabel(), pos, sequence);
                            System.out.format("%s\t%s\t%s%n", sensor.getMeta(), sensor.getCodon(),
                                    sensor.toString());
                        }
                    }
                }
                inStream.close();
            }
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }

    }

}
