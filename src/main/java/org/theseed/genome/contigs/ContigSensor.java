/**
 *
 */
package org.theseed.genome.contigs;

import java.util.Arrays;
import java.util.List;
import java.util.stream.Collectors;

/**
 * A contig sensor represents a position on a contig, and contains information about
 * the DNA letters at and to either side of that position.  Each DNA letter is represented
 * by one or more floating-point numbers.  A ContigSensorFactory class builds these.
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

    /**
     * Construct a contig sensor for a specified position on a contig.
     *
     * @param id	ID of the contig
     * @param pos	position of the sensor
     */
    public ContigSensor(String id, int pos) {
        this.contigId = id;
        this.position = pos;
        this.suspicious = false;
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
     * @return the sensor values
     */
    public double[] getSensors() {
        return this.sensors;
    }

    /**
     * Store the sensor values in this object.
     *
     * @param newSensors	new sensors to store
     * @param suspicious	TRUE if there were ambiguity characters found
     */
    protected void storeSensors(double[] newSensors, boolean suspicious) {
        this.sensors = newSensors;
        this.suspicious = suspicious;
    }

    /**
     * @return the sensor values in a list
     */
    public List<Double> getSensorList() {
        List<Double> retVal = Arrays.stream(this.sensors).boxed().collect(Collectors.toList());
        return retVal;
    }

    /**
     * @return TRUE if this sensor includes ambiguity characters
     */
    public boolean isSuspicious() {
        return this.suspicious;
    }

    /**
     * @return the number of sensor columns in this object
     */
    public Object getSensorWidth() {
        return this.sensors.length;
    }


}
