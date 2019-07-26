/**
 *
 */
package org.theseed.genome.contigs;

/**
 * This factory assigns a floating-point number to each DNA letter and presents each letter as a single
 * input column.  The coding is A = -0.3, C = -0.6, G = 0.6, and T = 0.3.  Ambiguity characters are 0.
 * @author Bruce Parrello
 *
 */
public class DirectContigSensorFactory extends ContigSensorFactory {

    @Override
    protected void convertSequence(ContigSensor sensor, String sequence, int pos) {
        boolean suspicion = false;
        int halfWidth = ContigSensorFactory.getHalfWidth();
        int fullWidth = halfWidth * 2 + 1;
        int offset = pos - halfWidth - 1;
        double[] buffer = new double[fullWidth];
        for (int i = 0; i < buffer.length; i++) {
            int actual = offset + i;
            if (actual < 0 || actual >= sequence.length()) {
                buffer[i] = 0.0;
            } else switch (sequence.charAt(actual)) {
                case 'A' :
                case 'a' :
                    buffer[i] = -0.3;
                    break;
                case 'C' :
                case 'c' :
                    buffer[i] = -0.6;
                    break;
                case 'G' :
                case 'g' :
                    buffer[i] = 0.6;
                    break;
                case 'T' :
                case 't' :
                    buffer[i] = 0.3;
                    break;
                default :
                    buffer[i] = 0.0;
                    suspicion = true;
            }
        }
        sensor.storeSensors(buffer, suspicion);
    }

}
