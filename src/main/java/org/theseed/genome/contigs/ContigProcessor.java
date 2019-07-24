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
import org.theseed.counters.CountMap;
import org.theseed.genomes.Contig;
import org.theseed.genomes.Genome;
import org.theseed.genomes.GenomeDirectory;
import org.theseed.locations.Frame;
import org.theseed.locations.LocationList;
import org.theseed.utils.ICommand;

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
 * 		the contig this size or greater; the default is 90000
 * -f	forward only; the sensors are only to the right of the target position
 * -n	normally, only plus-strand locations are considered coding regions; if this is
 * 		specified, minus-strand locations are included as well
 *
 * --type		type of classification to do; the values are
 *    	coding	outputs a class of "coding" for a frame in a coding region and
 * 				"space" for a frame not in a coding region; the default is to
 * 				output the actual frame label
 * 		edge	outputs a class of "start" for the first base pair of a coding region,
 * 				"stop" for the last base pair of a coding region, and "space" for
 * 				everything else
 * 		phase	outputs a class of "0" for a non-coding region, "+1" for the first
 * 				base pair of a codon, "+2" for the second, and "+3" for the third
 * 				(this is the default)
 *
 * The positional parameters are the names of the input directories.
 *
 * @author Bruce Parrello
 *
 */

public class ContigProcessor implements ICommand {

    // FIELDS
    /** random number generator */
    private static Random rand = new Random();
    /** tracker for the number of examples generated per frame */
    private CountMap<String> classCounter;

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

    /** forward-only flag */
    @Option(name="-f", aliases={"--rightOnly", "--forwardOnly"},
            usage="if specified, no sensors left of the position are used")
    private void setPlusOnly(boolean newFlag) {
        ContigSensor.setPlusOnly(newFlag);
    }

    /** run length for output */
    @Option(name="-r", aliases={"--run"}, metaVar="200", usage="length of an output run per contig")
    private int runLength;

    /** debug switch */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="write progress messages to STDERR")
    private boolean debug;

    /** contig chunk size */
    @Option(name="-k", aliases={"--chunk"}, metaVar="50000",
            usage="size of a contig section for choosing a run")
    private int chunkSize;

    /** negative-allowed flag */
    @Option(name="-n", aliases= {"--negative", "--minus"}, usage="include minus strand results")
    private boolean negative;

    /** output type flag */
    @Option(name="--type", usage="type of classification")
    private LocationClass.Type classType;

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
        this.help = false;
        this.debug = false;
        this.chunkSize = 90000;
        this.negative = false;
        this.classType = LocationClass.Type.PHASE;
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
        // Initialize the private data.
        this.classCounter = new CountMap<String>();
        LocationClass lsensor = LocationClass.scheme(this.classType, this.negative);
        // The first job is to create the output header.  The first column is the
        // frame and the remaining columns are sensors.
        System.out.print("frame");
        ContigSensor.sensor_headers(System.out);
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
                        ProcessContig(contig, codingMap.get(contig.getId()), lsensor);
                    }
                }
            }
            if (debug) {
                // Display counts for each frame, so we can see if we have well-distributed
                // results.
                for (String cl : this.classCounter.keys()) {
                    System.err.format("%s has %8d results.%n", cl,
                            this.classCounter.getCount(cl));
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
     * @param lsensor 	classification scheme for locations
     */
    private void ProcessContig(Contig contig, LocationList framer, LocationClass lsensor) {
        // Run through the contig in 50,000 base-pair chunks, choosing random locations
        // to output.
        int limit = contig.length();
        int pos = 1;
        int end = pos + this.chunkSize;
        while (pos <= limit) {
            // Insure the end is not too big.
            if (end > limit) end = limit + 1;
            // Choose a random number >= the position and < end.
            int start = rand.nextInt(end - pos) + pos;
            // This will count the number of valid positions output.
            int count = 0;
            // This tracks the previous frame.
            Frame prevFrame = Frame.P0;
            // Loop through the contig locations.
            while (start <= limit && count < this.runLength) {
                // Get the coding frame. That's our first output column.
                Frame thisFrame = framer.computeRegionFrame(start, start);
                // Try to output this position.  We only output non-suspicious positions for
                // valid frames.
                if (thisFrame != Frame.XX) {
                    ContigSensor proposal = new ContigSensor(contig.getId(), start, contig.getSequence());
                    if (! proposal.isSuspicious()) {
                        // Compute the frame string.
                        String frame = lsensor.classOf(prevFrame, thisFrame);
                        // Write the frame followed by the sensor data.
                        System.out.format("%s\t%s%n", frame, proposal.toString());
                        // Record the output.
                        count++;
                        this.classCounter.count(frame);
                    }
                }
                // Move to the next position.
                start++;
                prevFrame = thisFrame;
            }
            pos = (start >= end ? start + 1 : end);
            end = pos + this.chunkSize;
        }
    }

}
