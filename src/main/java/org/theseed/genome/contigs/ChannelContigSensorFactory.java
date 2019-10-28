/**
 *
 */
package org.theseed.genome.contigs;

/**
 * This very simple sensor factory simply copies each DNA letter to the output as a column.
 * We have a rather extensive switch statement to insure only valid letters are copied.
 *
 * @author Bruce Parrello
 *
 */
public class ChannelContigSensorFactory extends ContigSensorFactory {

    @Override
    protected void convertSequence(ContigSensor sensor, String sequence, int pos) {
        boolean suspicion = false;
        int offset = pos - ContigSensorFactory.getLeftWidth() - 1;
        String[] buffer = new String[ContigSensorFactory.getFullWidth()];
        for (int i = 0; i < buffer.length; i++) {
            int actual = offset + i;
            if (actual < 0 || actual >= sequence.length()) {
                buffer[i] = "-";
            } else {
                char nucleon = sequence.charAt(actual);
                switch (nucleon) {
                case 'a' :
                case 'A' :
                case 'c' :
                case 'C' :
                case 'g' :
                case 'G' :
                case 't' :
                case 'T' :
                case 'u' :
                case 'U' :
                case '-' :
                case 'x' :
                case 'X' :
                case 'y' :
                case 'Y' :
                case 'r' :
                case 'R' :
                case 'w' :
                case 'W' :
                case 's' :
                case 'S' :
                case 'k' :
                case 'K' :
                case 'm' :
                case 'M' :
                    buffer[i] = String.valueOf(nucleon);
                    break;
                default :
                    buffer[i] = "X";
                    suspicion = true;
                }
            }
        }
        sensor.storeSensors(buffer, suspicion);
    }

}
