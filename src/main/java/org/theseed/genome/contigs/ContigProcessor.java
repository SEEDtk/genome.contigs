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
import org.theseed.io.BalancedOutputStream;
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
 * -n	normally, only plus-strand locations are considered coding regions; if this is
 * 		specified, minus-strand locations are included as well
 * -f	filter on known edge codons
 * -b	indicates the output should be balanced; the records are held in memory and then
 * 		a distributed, balanced subset is output; the value should be a number from 1.0
 * 		to 2.0, indicating the maximum number of output records per class as a fraction of
 * 		the smallest class's size
 *
 * --type		type of classification to do; the values are
 *    	coding	outputs a class of "coding" for a frame in a coding region and
 * 				"space" for a frame not in a coding region; the default is to
 * 				output the actual frame label
 * 		edge	outputs a class of "start" for the first base pair of a coding region,
 * 				"stop" for the last base pair of a coding region, and "other" for
 * 				everything else (this is the default)
 * 		start	outputs a class of "start" for the first base pair of a coding region,
 * 				and "other" for everything else
 * 		stop	outputs a class of "stop" for the first base pair past the end of a coding region,
 * 				and "other" for everything else
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
    /** factory object for creating contig sensors */
    private ContigSensorFactory factory;
    /** data output stream */
    private BalancedOutputStream outStream;

    // COMMAND-LINE OPTIONS

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** sensor width */
    @Option(name="-w", aliases={"--width"}, metaVar="14",
            usage="distance on either side for sensors")
    private void setWidth(int newWidth) {
        ContigSensorFactory.setHalfWidth(newWidth);
    }

    /** balanced output fuzz factor */
    @Option(name="-b", aliases={"--balance", "--fuzz"}, metaVar="1.2", usage="specify class-balanced output")
    private double fuzzFactor;

    /** filter for edge codons */
    @Option(name="-f", aliases={"--edgeFilter"}, usage="codon filtering type")
    private boolean edgeFilter;

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

    /** sensor type */
    @Option(name="--sensor", metaVar="one_hot", usage="type of DNA sensor to use")
    private void setFactory(ContigSensorFactory.Type type) {
        this.factory = ContigSensorFactory.create(type);
    }

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
        this.classType = LocationClass.Type.EDGE;
        this.edgeFilter = false;
        this.fuzzFactor = 0;
        this.factory = ContigSensorFactory.create(ContigSensorFactory.Type.CHANNEL);
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
                // Validate the fuzz factor.
                if (this.fuzzFactor != 0 && (this.fuzzFactor < 1.0 || this.fuzzFactor > 2.0)) {
                    throw new IllegalArgumentException("Balance factor must be 0 (off) or between 1.0 and 2.0 inclusive.");
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
        // Create the output stream.
        this.outStream = new BalancedOutputStream(this.fuzzFactor, System.out);
        // Set up the edge filter.
        CodonFilter filter = null;
        if (this.edgeFilter)
            filter = LocationClass.filter(this.classType);
        // The first job is to create the output header.  The first column is the
        // frame and the remaining columns are sensors.
        this.outStream.writeImmediate("frame", this.factory.sensor_headers());
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
                        ProcessContig(contig, codingMap.get(contig.getId()), lsensor, filter);
                    }
                }
            }
            this.outStream.close();
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
        // Close the output stream.  This is where the IO happens.
    }

    /**
     * Output the training data from the specified contig.
     *
     * @param contig	contig of interest
     * @param framer	location list used to compute frames
     * @param lsensor 	classification scheme for locations
     * @param filter	optional codon filter
     */
    private void ProcessContig(Contig contig, LocationList framer, LocationClass lsensor, CodonFilter filter) {
        // Activate the contig.
        lsensor.setLocs(framer);
        // Run through the contig in chunks, choosing random locations to output.
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
            // Extract the contig sequence.
            String sequence = contig.getSequence();
            // Loop through the contig locations.
            while (start <= limit && count < this.runLength) {
                if (filter == null || filter.matches(start, sequence)) {
                    ContigSensor proposal = this.factory.create(contig.getId(), start, sequence);
                    if (! proposal.isSuspicious()) {
                        // Compute the frame string.
                        String frame = lsensor.classOf(start);
                        if (frame != null) {
                            // Write the frame followed by the sensor data.
                            this.outStream.write(frame, proposal.toString());
                            // Record the output.
                            count++;
                            this.classCounter.count(frame);
                        }
                    }
                }
                // Move to the next position.
                start++;
            }
            pos = (start >= end ? start + 1 : end);
            end = pos + this.chunkSize;
        }
    }

}
