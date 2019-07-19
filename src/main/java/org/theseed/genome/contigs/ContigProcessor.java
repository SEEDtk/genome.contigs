/**
 *
 */
package org.theseed.genome.contigs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.List;
import java.util.Map;
import java.util.Random;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genomes.Contig;
import org.theseed.genomes.Genome;
import org.theseed.genomes.GenomeDirectory;
import org.theseed.locations.LocationList;

/**
 * This procedure reads all the genomes in a GTO directory and produces a training
 * set for computing coding frames.  In each contig processed, we choose a random
 * region and output a run of positions along with their coding frames.  For each
 * position, the nucleotides for a certain distance on either side (the "half-width")
 * are encoded into floating-point numbers.  The output file will be tab-delimited
 * with headers.
 *
 * The following command-line options are supported.
 *
 * -w	the number of positions to look at on either side of a base pair; the default
 * 		is 14
 * -r	the number of positions in a row to output; the default is 200
 * -v	write progress messages to STDERR
 * -k	chunk size for contigs; one output run will be produced for every piece of
 * 		the contig this size or greater; the default is 50000
 *
 * The positional parameters are the names of the input directories.
 *
 * @author Bruce Parrello
 *
 */
public class ContigProcessor {

    // COMMAND-LINE OPTIONS

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** sensor width */
    @Option(name="-w", aliases={"--width"}, metaVar="14",
            usage="distance on either side for sensors")
    private void setWidth(int newWidth) {
        ContigSensor.setHalfWidth(newWidth);
    }

    /** run length for output */
    @Option(name="-r", aliases={"--run"}, metaVar="200", usage="length of an output run per contig")
    private int runLength;

    /** debug switch */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="write progress messages to STDERR")
    private boolean debug;

    /** contig chunk size */
    @Option(name="-k", aliases={"--chunk"}, usage="size of a contig section for choosing a run")
    private int chunkSize;

    /** input directories */
    @Argument(index=0, metaVar="genomeDir1 genomeDir2 ...", usage="directories containing genome objects",
            multiValued=true)
    private List<File> genomeDirs;

    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.runLength = 200;
        this.debug = false;
        this.chunkSize = 50000;
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // Insure the genome directories are valid.
                for (File genomeDir : genomeDirs) {
                    if (! genomeDir.isDirectory()) {
                        throw new FileNotFoundException(genomeDir.getPath() + " is not a valid directory.");
                    }
                }
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println("Error processing genome directory: " + e.getMessage());
        }
        return retVal;
    }


    /**
     * Process the genome directories to produce the output file.
     */
    public void run() {
        // The first job is to create the output header.  The first column is the
        // frame and the remaining columns are sensors.  We label the sensors "pos."
        // and then the number "-width" to "+width".
        System.out.print("frame");
        int hw = ContigSensor.getHalfWidth();
        for (int idx = -hw; idx <= hw; idx++) {
            System.out.format("\tpos.%02d", idx);
        }
        System.out.println();
        try {
            // Loop through the genome directories.
            for (File genomeDir : this.genomeDirs) {
                if (debug) System.err.println("Processing " + genomeDir + ".");
                GenomeDirectory genomes = new GenomeDirectory(genomeDir.getPath());
                // Loop through the genomes.
                for (Genome genome : genomes) {
                    if (debug) System.err.println("Processing " + genome + ".");
                    // Create this genome's coding map.
                    Map<String, LocationList> codingMap = LocationList.createGenomeCodingMap(genome);
                    for (Contig contig : genome.getContigs()) {
                        ProcessContig(contig, codingMap.get(contig.getId()));
                    }
                }
            }
        } catch (IOException e) {
            System.err.println("Error processing genome directory: " + e.getMessage());
        }
    }


    /**
     * Output the training data from the specified contig.
     *
     * @param contig	contig of interest
     * @param framer	location list used to compute frames
     */
    private void ProcessContig(Contig contig, LocationList framer) {
        // Get a random number generator.
        Random rand = new Random();
        // Run through the contig in 50,000 base-pair chunks, choosing random locations
        // to output.
        int limit = contig.length();
        int pos = 1;
        int end = pos + this.chunkSize;
        while (pos <= limit) {
            // TODO pick and output a run
            pos = end;
            end = pos + this.chunkSize;
        }

    }

}
