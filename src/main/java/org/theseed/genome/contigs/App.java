package org.theseed.genome.contigs;

import java.util.Arrays;

import org.theseed.utils.ICommand;

/**
 * Output contig data for the learning module.  The possible commands are "train" to output
 * a training set, "predict" to output an input set for prediction, and "test" to
 * output a test set for verification.
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        ICommand processor;
        // Parse the parameters.
        switch (command) {
        case "train" :
            processor = new ContigProcessor();
            break;
        case "predict" :
            processor = new FastaProcessor();
            break;
        case "test" :
            processor = new GenomeProcessor();
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ": must be \"train\", \"test\", or \"predict\".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }

}
