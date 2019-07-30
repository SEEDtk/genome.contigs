/**
 *
 */
package org.theseed.genome.contigs;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.counters.CountMap;
import org.theseed.genomes.Contig;
import org.theseed.genomes.Genome;
import org.theseed.locations.Frame;
import org.theseed.locations.LocationList;
import org.theseed.utils.ICommand;

/**
 * This procedure reads a single genome GTO file and generates a verification
 * file for the deep learning engine.  Every base pair is converted into a
 * row of the output file, along with metadata indicating its location and
 * the expected value.  When the file is input to the model, its predictions
 * can be compared to the expected value in the metadata.
 *
 * -v	write progress messages to STDERR
 * -n	normally, only plus-strand locations are considered coding regions; if this is
 * 		specified, minus-strand locations are included as well
 *
 * --type		type of classification to do; the values are
 *    	coding	outputs a class of "coding" for a frame in a coding region and
 * 				"space" for a frame not in a coding region; the default is to
 * 				output the actual frame label
 * 		edge	outputs a class of "start" for the first base pair of a coding region,
 * 				"stop" for the last base pair of a coding region, and "space" for
 * 				everything else (this is the default)
 * 		phase	outputs a class of "0" for a non-coding region, "+1" for the first
 * 				base pair of a codon, "+2" for the second, and "+3" for the third
 *
 * --sensor		type of DNA sensor to use
 * 		direct	each base pair converts to a single number
 * 		codon	each base pair converts to a number computed from the three base pairs beginning
 * 				at the current position
 * 		channel	each base pair is converted to a vector of the probabilities for each possible
 * 				nucleotide value (this is the default)
 *
 * This positional parameter is the name of the GTO file containing the genome.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProcessor implements ICommand {

    // FIELDS
    /** factory object for creating contig sensors */
    private ContigSensorFactory factory;

    // COMMAND-LINE OPTIONS

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** debug switch */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="write progress messages to STDERR")
    private boolean debug;

    /** negative-allowed flag */
    @Option(name="-n", aliases= {"--negative", "--minus"}, usage="include minus strand results")
    private boolean negative;

    /** output type flag */
    @Option(name="--type", usage="type of classification")
    private LocationClass.Type classType;

    /** sensor type */
    @Option(name="--sensor", metaVar="codon", usage="type of DNA sensor to use (default: CHANNEL)")
    private void setFactory(ContigSensorFactory.Type type) {
        this.factory = ContigSensorFactory.create(type);
    }

    /** input file */
    @Argument(index=0, metaVar="genomeFile", usage="file containing the genome object")
    private File genomeFile;

    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.debug = false;
        this.negative = false;
        this.classType = LocationClass.Type.EDGE;
        this.factory = ContigSensorFactory.create(ContigSensorFactory.Type.CHANNEL);
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // Insure the genome file is valid.
                if (! this.genomeFile.exists()) {
                    throw new FileNotFoundException(genomeFile + " does not exist.");
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
        // Create the location sensor.
        LocationClass lsensor = LocationClass.scheme(this.classType, this.negative);
        // Read in the genome.
        try {
            Genome genome = new Genome(genomeFile);
            // Create the output header.  The first column is the
            // location, then the expection,  and finally the sensors.
            System.out.println("location\texpect\t" + this.factory.sensor_headers());
            // We use this to count the output classes.
            CountMap<String> classCounts = new CountMap<String>();
            // Get the genome's contig map.
            Map<String, LocationList> codingMap = LocationList.createGenomeCodingMap(genome);
            for (Contig contig : genome.getContigs()) {
                if (debug) System.err.println("Processing contig " + contig.getId());
                // Get this contig's location list.
                LocationList framer = codingMap.get(contig.getId());
                // Loop through the base pairs, generating data.
                int limit = contig.length();
                // This tracks the previous frame.
                Frame prevFrame = Frame.P0;
                // Loop through the contig.
                for (int pos = 1; pos <= limit; pos++) {
                    Frame thisFrame = framer.computeRegionFrame(pos, pos);
                    // Compute this location's expected value. If the expectation
                    // is unknown, we use a question mark.
                    String expect = lsensor.classOf(prevFrame, thisFrame);
                    if (expect == null) expect = "?";
                    classCounts.count(expect);
                    // Compute this location's sensor values.
                    ContigSensor proposal = this.factory.create(contig.getId(), pos, contig.getSequence());
                    // Write it all out.
                    System.out.format("%s\t%s\t%s%n", proposal.getMeta(),
                            expect, proposal.toString());
                    // Set up for the next iteration.
                    prevFrame = thisFrame;
                }
            }
            if (this.debug) for (CountMap<String>.Count count : classCounts.sortedCounts()) {
                System.err.format("%20d written of type %s%n", count.getCount(), count.getKey());
            }
   } catch (NumberFormatException | IOException e) {
            System.err.println("Error processing " + genomeFile + ": " +
                    e.getMessage());
        }
    }

}
