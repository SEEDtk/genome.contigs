/**
 *
 */
package org.theseed.genome.contigs;

import java.util.ArrayList;

import org.apache.commons.lang3.StringUtils;

/**
 * This factory assigns four input channels to each DNA position.  Three of them will have a value of 0
 * and the fourth will be 1.  The channel chosen depends on the nucleotide.
 *
 * @author Bruce Parrello
 *
 */
public class OneHotContigSensorFactory extends ContigSensorFactory {

    /** source arrays with the channel values */
    public static final double[] A_VALUES = { 1.0, 0.0, 0.0, 0.0 };
    public static final double[] C_VALUES = { 0.0, 1.0, 0.0, 0.0 };
    public static final double[] G_VALUES = { 0.0, 0.0, 1.0, 0.0 };
    public static final double[] T_VALUES = { 0.0, 0.0, 0.0, 1.0 };
    public static final double[] X_VALUES = { 0.0, 0.0, 0.0, 0.0 };

    @Override
    public String sensor_headers() {
        int halfWidth = ContigSensorFactory.getHalfWidth();
        int fullWidth = halfWidth * 2 + 1;
        ArrayList<String> headers = new ArrayList<String>(fullWidth);
        for (int i = -halfWidth; i <= halfWidth; i++) {
            headers.add(String.format("pos.%dA\tpos.%dC\tpos.%dG\tpos.%dT", i, i, i, i));
        }
        return StringUtils.join(headers, '\t');
    }

    @Override
    protected void convertSequence(ContigSensor sensor, String sequence, int pos) {
        boolean suspicion = false;
        int halfWidth = ContigSensorFactory.getHalfWidth();
        int fullWidth = halfWidth * 2 + 1;
        double[] buffer = new double[fullWidth * 4];
        int offset = pos - halfWidth - 1;
        int idx = 0;
        for (int i = 0; i < fullWidth; i++) {
            int actual = offset + i;
            double[] chosen = X_VALUES;
            if (actual >= 0 && actual < sequence.length()) {
                switch (sequence.charAt(actual)) {
                case 'A' :
                case 'a' :
                    chosen = A_VALUES;
                    break;
                case 'C' :
                case 'c' :
                    chosen = C_VALUES;
                    break;
                case 'G' :
                case 'g' :
                    chosen = G_VALUES;
                    break;
                case 'T' :
                case 't' :
                case 'U' :
                case 'u' :
                    chosen = T_VALUES;
                    break;
                default :
                    suspicion = true;
                }
            }
            for (double v : chosen) {
                buffer[idx++] = v;
            }
        }
        sensor.storeSensors(buffer, suspicion);
    }

}
