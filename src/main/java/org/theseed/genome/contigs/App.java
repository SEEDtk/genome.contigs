package org.theseed.genome.contigs;

import java.util.Arrays;

/**
 * Output contig data for the learning module.  The possible commands are "train" to output
 * a training set and "predict" to output an input set for prediction.  The training set
 * output is handled by ContigProcessor and the prediction set output by SequenceProcessor.
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        // Parse the parameters.
        switch (command) {
        case "train" :
            ProcessContigs(newArgs);
            break;
        default :
            System.err.println("Invalid command " + command + ".");
        }
    }

    /**
     * Process contigs to produce training output.
     *
     * @param newArgs	command-line parameters
     */
    private static void ProcessContigs(String[] newArgs) {
        ContigProcessor processor = new ContigProcessor();
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }

}
