/**
 *
 */
package org.theseed.genome.contigs;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;
import org.theseed.genomes.Contig;
import org.theseed.sequence.Sequence;

/**
 * A contig sensor represents a position on a contig, and contains information about
 * the DNA letters at and to either side of that position.  In normal mode, the position
 * itself is in the exact center of the sensor.  In plus-only mode, the position is on the
 * left edge of the sensor.  Each DNA letter is represented by a floating-
 * point number, with other/missing = 0.
 *
 * Note that in Deep Learning lingo, what we are calling a "sensor" is called a
 * "feature".  That term, however, is already heavily overloaded on the genome
 * side of things.
 *
 * @author Bruce Parrello
 *
 */
public class ContigSensor {

    // FIELDS
    /** ID of the source contig */
    String contigId;
    /** position in the contig (1-based) */
    int position;
    /** sensors to either side */
    double[] sensors;
    /** TRUE if this sensor includes ambiguity characters */
    boolean suspicious;

    /** global sensor width, to either side of the target position */
    private static int halfWidth = 14;

    /** TRUE if we should only sense to the right */
    private static boolean plusOnly = false;

    /**
     * Construct a contig sensor from a specific position on the contig.
     *
     * @param contigId	ID of the contig
     * @param position	position (1-based) of the point of interest
     * @param sequence	DNA sequence of the contig
     */
    public ContigSensor(String contigId, int position, String sequence) {
        this.contigId = contigId;
        this.position = position;
        this.suspicious = false;
        // Now fill in the sensors.
        this.sensors = new double[getSensorWidth()];
        // Compute the first contig offset to be sensed.
        int offset0 = position - 1 - (plusOnly ? 0 : halfWidth);
        for (int i = 0; i < this.sensors.length; i++) {
            int p = offset0 + i;
            if (p < 0 || p >= sequence.length()) {
                this.sensors[i] = 0.0;
            } else {
                switch (sequence.charAt(p)) {
                case 'A':
                case 'a' :
                    this.sensors[i] = -0.3;
                    break;
                case 'C':
                case 'c':
                    this.sensors[i] = -0.6;
                    break;
                case 'G':
                case 'g':
                    this.sensors[i] = 0.6;
                    break;
                case 'T':
                case 't':
                case 'U':
                case 'u':
                    this.sensors[i] = 0.3;
                    break;
                default:
                    this.sensors[i] = 0.0;
                    this.suspicious = true;
                }
            }
        }
    }

    /**
     * @return the contigId for this sensor
     */
    public String getContigId() {
        return this.contigId;
    }

    /**
     * @return a tab-delimited string of all the sensor values, in order
     */
    @Override
    public String toString() {
        return Arrays.stream(this.sensors).mapToObj(Double::toString)
                .collect(Collectors.joining("\t"));
    }
    /**
     * @return the position for this sensor
     */
    public int getPosition() {
        return this.position;
    }

    /**
     * @return a string representation of contig and position
     */
    public String getMeta() {
        return this.contigId + ";" + Integer.toString(this.position);
    }

    /**
     * @return the sensors
     */
    public double[] getSensors() {
        return this.sensors;
    }

    /**
     * @return the sensors in a list
     */
    public List<Double> getSensorList() {
        ArrayList<Double> retVal = new ArrayList<Double>(getSensorWidth());
        for (double item : this.sensors) {
            retVal.add(item);
        }
        return retVal;
    }

    /**
     * @return TRUE if this sensor includes ambiguity characters
     */
    public boolean isSuspicious() {
        return this.suspicious;
    }

    /**
     * @return the number of sensors on each side of the target position
     */
    public static int getHalfWidth() {
        return halfWidth;
    }

    /**
     * Specify a new global half-width.  This does not affect sensors
     * already constructed.
     *
     * @param new sensor width
     */
    public static void setHalfWidth(int newWidth) {
        ContigSensor.halfWidth = newWidth;
    }

    /**
     * @return the number of values in a sensor.
     */
    public static int getSensorWidth() {
        return halfWidth + 1 + (plusOnly ? 0 : halfWidth);
    }

    /**
     * @return all the non-suspicious sensors for the specified FASTA sequence
     *
     * @param contigSequence	FASTA sequence representing a contig
     */
    public static List<ContigSensor> processContig(Sequence contigSequence) {
        return processContig(contigSequence.getLabel(), contigSequence.getSequence(),
                1, contigSequence.length());
    }

    /**
     * @return all the non-suspicious sensors for the specified Contig
     *
     * @param contig	contig containing the sequence to process
     * @param start		starting position on the contig (1-based)
     * @param len		number of positions on the contig to return
     */
    public static List<ContigSensor> processContig(Contig contig, int start, int len) {
        return processContig(contig.getId(), contig.getSequence(), start, len);
    }

    /**
     * @return all the non-suspicious sensors for an identified contig
     *
     * @param contigId	ID of the contig
     * @param sequence	sequence of the contig
     * @param start		starting position on the contig (1-based)
     * @param len		number of positions on the contig to return
     */
    private static List<ContigSensor> processContig(String contigId, String sequence,
            int start, int len) {
        ArrayList<ContigSensor> retVal = new ArrayList<ContigSensor>(sequence.length());
        int end = start + len - 1;
        if (end > sequence.length()) end = sequence.length();
        for (int i = start; i <= end; i++) {
            ContigSensor snapshot = new ContigSensor(contigId, i, sequence);
            if (! snapshot.suspicious) retVal.add(snapshot);
        }
        return retVal;
    }

    /**
     * Write the sensor headers to the specified stream.  It is presumed they are not
     * in the first column.
     *
     * @param out	desired output stream
     */
    public static void sensor_headers(PrintStream out) {
        int hw = ContigSensor.getHalfWidth();
        for (int idx = (plusOnly ? 0 : -hw); idx <= hw; idx++) {
            out.format("\tpos.%02d", idx);
        }
    }

    /**
     * @return TRUE if the sensors are only to the right of the target position
     */
    public static boolean isPlusOnly() {
        return plusOnly;
    }

    /**
     * @param plusOnly 	TRUE if the sensors should only be to the right of the target position,
     * 					else FALSE
     */
    public static void setPlusOnly(boolean plusOnly) {
        ContigSensor.plusOnly = plusOnly;
    }



}
