/**
 *
 */
package org.theseed.genome.contigs;

/**
 * In this class, the sensor value is computed from the three nucleotides at the current position.
 * Each nucleotide corresponds to a digit past the decimal point, and the values are A = 0.2, C = 0.6,
 * G = 0.8, and T = 0.4. So (for example) ACG would be 0.268.
 *
 * @author Bruce Parrello
 *
 */
public class CodonContigSensorFactory extends ContigSensorFactory {

    @Override
    protected void convertSequence(ContigSensor sensor, String sequence, int pos) {
        boolean suspicion = false;
        int halfWidth = ContigSensorFactory.getHalfWidth();
        int fullWidth = halfWidth * 2 + 1;
        int offset = pos - halfWidth - 1;
        String[] buffer = new String[fullWidth];
        for (int i = 0; i < buffer.length; i++) {
            int actual = offset + i;
            StringBuffer v = new StringBuffer(5);
            v.append("0.");
            for (int j = actual; j <= actual + 2; j++) {
                if (j >= 0 && j < sequence.length()) {
                    switch (sequence.charAt(j)) {
                    case 'A' :
                    case 'a' :
                        v.append('2');
                        break;
                    case 'C' :
                    case 'c' :
                        v.append('6');
                        break;
                    case 'G' :
                    case 'g' :
                        v.append('8');
                        break;
                    case 'T' :
                    case 't' :
                    case 'U' :
                    case 'u' :
                        v.append('4');
                        break;
                    default :
                        v.append('0');
                        suspicion = true;
                    }
                }
            }
            buffer[i] = v.toString();
        }
        sensor.storeSensors(buffer, suspicion);
    }

}
