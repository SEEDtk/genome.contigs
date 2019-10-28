package org.theseed.genome.contigs;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.util.List;

import org.theseed.genome.Contig;
import org.theseed.locations.Frame;
import org.theseed.locations.Location;
import org.theseed.locations.LocationList;
import org.theseed.sequence.Sequence;



/**
 * Unit test for simple App.
 */
public class AppTest
    extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }


    /**
     * Test the contig sensors.
     */
    public void testContigSensors()
    {
        String contigID = "3000.contig.1";
        Sequence frec = new Sequence(contigID, "", "AACGTCCTGAAGTC");
        ContigSensorFactory.setLeftWidth(4);
        ContigSensorFactory.setRightWidth(4);
        ContigSensorFactory myFactory = ContigSensorFactory.create(ContigSensorFactory.Type.DIRECT);
        assertThat("Wrong header", myFactory.sensor_headers(),
                equalTo("pos.-4\tpos.-3\tpos.-2\tpos.-1\tpos.0\tpos.1\tpos.2\tpos.3\tpos.4"));
        List<ContigSensor> sensors = myFactory.processContig(frec);
        ContigSensor sensor1 = sensors.get(0);
        ContigSensor sensor11 = sensors.get(11);
        assertThat("Wrong contig ID for first sensor.", sensor1.getContigId(), equalTo(contigID));
        assertThat("Wrong position for first sensor.", sensor1.getPosition(), equalTo(1));
        assertThat("Wrong meta string for first sensor.", sensor1.getMeta(), equalTo(contigID + ";1"));
        assertThat("Wrong sensors for first", sensor1.getSensorList(),
                contains("0.0", "0.0", "0.0", "0.0", "-0.3", "-0.3", "-0.6", "0.6", "0.3"));
        assertThat("Wrong string", sensor1.toString(), equalTo("0.0\t0.0\t0.0\t0.0\t-0.3\t-0.3\t-0.6\t0.6\t0.3"));
        assertThat("Wrong sensors for 11th", sensor11.getSensorList(),
                contains("0.3", "0.6", "-0.3", "-0.3", "0.6", "0.3", "-0.6", "0.0", "0.0"));
        assertThat("Wrong contig ID for 11th", sensor11.getContigId(), equalTo(contigID));
        assertThat("Wrong position for 11th", sensor11.getPosition(), equalTo(12));
        for (ContigSensor sensor : sensors) {
            assertThat("Wrong contig ID in sensor", sensor.getContigId(), equalTo(contigID));
        }
        Contig contig = new Contig(contigID, "AACGTNCGGGGAAAT", 11);
        sensors = myFactory.processContig(contig, 2, 20);
        sensor11 = sensors.get(0);
        assertThat("Wrong position for post-N sensor.", sensor11.getPosition(), equalTo(11));
        sensor1 = myFactory.create(contigID, 7, contig.getSequence());
        assertTrue("Sensor not suspicious.", sensor1.isSuspicious());
        String[] sensorArray = sensor1.getSensors();
        List<String> sensorList = sensor1.getSensorList();
        for (int i = 0; i < sensorArray.length; i++) {
            assertThat("Sensor mismatch at " + i, sensorArray[i], equalTo(sensorList.get(i)));
        }
        myFactory = ContigSensorFactory.create(ContigSensorFactory.Type.CHANNEL);
        sensors = myFactory.processContig(frec);
        sensor1 = sensors.get(0);
        sensor11 = sensors.get(10);
        assertThat("Wrong contig ID for first sensor.", sensor1.getContigId(), equalTo(contigID));
        assertThat("Wrong position for first sensor.", sensor1.getPosition(), equalTo(1));
        assertThat("Wrong meta string for first sensor.", sensor1.getMeta(), equalTo(contigID + ";1"));
        assertThat("Wrong sensors for first", sensor1.getSensorList(),
                contains("-", "-", "-", "-", "A", "A", "C", "G", "T"));
        assertThat("Wrong sensors for 11th", sensor11.getSensorList(),
                contains("C", "T", "G", "A", "A", "G", "T", "C", "-"));
        assertThat("Wrong contig ID for 10th", sensor11.getContigId(), equalTo(contigID));
        assertThat("Wrong position for 10th", sensor11.getPosition(), equalTo(11));
        myFactory = ContigSensorFactory.create(ContigSensorFactory.Type.CODON);
        ContigSensorFactory.setLeftWidth(3);
        ContigSensorFactory.setRightWidth(5);
        assertThat(ContigSensorFactory.getFullWidth(), equalTo(9));
        assertThat("Wrong header", myFactory.sensor_headers(),
                equalTo("pos.-3\tpos.0\tpos.3"));
        sensors = myFactory.processContig(frec);
        sensor1 = sensors.get(0);
        assertThat("Wrong contig ID for first sensor.", sensor1.getContigId(), equalTo(contigID));
        assertThat("Wrong position for first sensor.", sensor1.getPosition(), equalTo(1));
        assertThat("Wrong meta string for first sensor.", sensor1.getMeta(), equalTo(contigID + ";1"));
        assertThat("Wrong sensors for first", sensor1.getSensorList(),
                contains("---", "aac", "gtc"));
        assertThat("Wrong string", sensor1.toString(), equalTo("---\taac\tgtc"));
        sensor1 = sensors.get(1);
        assertThat("Wrong sensors for second", sensor1.getSensorList(),
                contains("--a", "acg", "tcc"));
        Sequence frec2 = new Sequence(contigID, "", "AACGTCCTRAAGTCA");
        myFactory = ContigSensorFactory.create(ContigSensorFactory.Type.AMINOACID);
        ContigSensorFactory.setLeftWidth(3);
        ContigSensorFactory.setRightWidth(5);
        assertThat(ContigSensorFactory.getFullWidth(), equalTo(9));
        assertThat("Wrong header", myFactory.sensor_headers(),
                equalTo("pos.-3\tpos.0\tpos.3"));
        sensors = myFactory.processContig(frec2);
        sensor1 = sensors.get(0);
        assertThat("Wrong contig ID for first sensor.", sensor1.getContigId(), equalTo(contigID));
        assertThat("Wrong position for first sensor.", sensor1.getPosition(), equalTo(1));
        assertThat("Wrong meta string for first sensor.", sensor1.getMeta(), equalTo(contigID + ";1"));
        assertThat("Wrong sensors for first", sensor1.getSensorList(),
                contains("-", "N", "V"));
        assertThat("Wrong string", sensor1.toString(), equalTo("-\tN\tV"));
        sensor1 = sensors.get(1);
        assertThat("Wrong sensors for second", sensor1.getSensorList(),
                contains("-", "T", "S"));
        sensor1 = sensors.get(3);
        assertThat("Failure of suspicion", sensor1.getSensorList(),
                contains("K", "S", "-"));

    }

    /*
     * Test location sensors.
     */
    public void testLocationClasses() {
        LocationClass lsensor = LocationClass.scheme(LocationClass.Type.CODING, true);
        assertThat("Incorrect normalize for M2, negative.", lsensor.normalize(Frame.M2), equalTo(Frame.M2));
        assertThat("Incorrect normalize for M1, negative.", lsensor.normalize(Frame.M1), equalTo(Frame.M1));
        assertThat("Incorrect normalize for M0, negative.", lsensor.normalize(Frame.M0), equalTo(Frame.M0));
        assertThat("Incorrect normalize for F0, negative.", lsensor.normalize(Frame.F0), equalTo(Frame.F0));
        assertThat("Incorrect normalize for P2, negative.", lsensor.normalize(Frame.P2), equalTo(Frame.P2));
        assertThat("Incorrect normalize for P1, negative.", lsensor.normalize(Frame.P1), equalTo(Frame.P1));
        assertThat("Incorrect normalize for P0, negative.", lsensor.normalize(Frame.P0), equalTo(Frame.P0));
        assertThat("Incorrect normalize for XX, negative.", lsensor.normalize(Frame.XX), equalTo(Frame.XX));
        LocationList newList = new LocationList("myContig");
        Location[] locs = { Location.create("myContig", "+", 10, 99, 102, 199),
                            Location.create("myContig", "-", 100, 400),
                            Location.create("myContig", "-", 500, 999),
                            Location.create("myContig", "-", 3000, 3199),
                            Location.create("myContig", "+", 4000, 4999),
                            Location.create("myContig", "-", 6000, 6199),
                            Location.create("myContig", "+", 6800, 7299),
                            Location.create("myContig", "+", 8000, 8199),
                            Location.create("myContig", "+", 8400, 8599),
                            Location.create("myContig", "+", 9100, 9190),
                            Location.create("myContig", "+", 9200, 9290),
                            Location.create("myContig", "-", 9300, 9390),
                            Location.create("myContig", "+", 9400, 9490),
                            Location.create("myContig", "-", 9500, 9999),
                          };
        // Add all these locations to the list.
        for (Location loc : locs) {
            assertTrue("Failed to add " + loc + " to list.", newList.addLocation(loc));
        }
        lsensor = LocationClass.scheme(LocationClass.Type.CODING, false);
        lsensor.setLocs(newList);
        assertNull(lsensor.classOf(150));
        assertThat(lsensor.classOf(4102), equalTo("coding"));
        assertThat(lsensor.classOf(6050), equalTo("space"));
        assertThat(lsensor.classOf(8000), equalTo("coding"));
        assertThat(lsensor.classOf(9999), equalTo("space"));
        assertThat(lsensor.classOf(9392), equalTo("space"));
        lsensor = LocationClass.scheme(LocationClass.Type.PHASE, true);
        lsensor.setLocs(newList);
        assertNull(lsensor.classOf(150));
        assertThat(lsensor.classOf(4102), equalTo("+1"));
        assertThat(lsensor.classOf(6050), equalTo("-3"));
        assertThat(lsensor.classOf(8000), equalTo("+1"));
        assertThat(lsensor.classOf(9999), equalTo("-1"));
        assertThat(lsensor.classOf(9392), equalTo("0"));
        lsensor = LocationClass.scheme(LocationClass.Type.PHASE, false);
        lsensor.setLocs(newList);
        assertNull(lsensor.classOf(150));
        assertThat(lsensor.classOf(4102), equalTo("+1"));
        assertThat(lsensor.classOf(6050), equalTo("0"));
        assertThat(lsensor.classOf(8000), equalTo("+1"));
        assertThat(lsensor.classOf(9999), equalTo("0"));
        assertThat(lsensor.classOf(9392), equalTo("0"));
        lsensor = LocationClass.scheme(LocationClass.Type.EDGE, false);
        lsensor.setLocs(newList);
        assertThat(lsensor.classOf(150), equalTo("other"));
        assertThat(lsensor.classOf(9999), equalTo("other"));
        assertThat(lsensor.classOf(8000), equalTo("start"));
        assertThat(lsensor.classOf(8597), equalTo("stop"));
        assertThat(lsensor.classOf(9499), equalTo("other"));
        lsensor = LocationClass.scheme(LocationClass.Type.EDGE, true);
        lsensor.setLocs(newList);
        assertThat(lsensor.classOf(150), equalTo("other"));
        assertThat(lsensor.classOf(9999), equalTo("start"));
        assertThat(lsensor.classOf(8000), equalTo("start"));
        assertThat(lsensor.classOf(8597), equalTo("stop"));
        assertThat(lsensor.classOf(9502), equalTo("stop"));
        lsensor = LocationClass.scheme(LocationClass.Type.START, false);
        lsensor.setLocs(newList);
        assertThat(lsensor.classOf(150), equalTo("other"));
        assertThat(lsensor.classOf(9999), equalTo("other"));
        assertThat(lsensor.classOf(8000), equalTo("start"));
        assertThat(lsensor.classOf(8597), equalTo("other"));
        assertThat(lsensor.classOf(9502), equalTo("other"));
        lsensor = LocationClass.scheme(LocationClass.Type.START, true);
        lsensor.setLocs(newList);
        assertThat(lsensor.classOf(150), equalTo("other"));
        assertThat(lsensor.classOf(9999), equalTo("start"));
        assertThat(lsensor.classOf(8000), equalTo("start"));
        assertThat(lsensor.classOf(8597), equalTo("other"));
        assertThat(lsensor.classOf(9502), equalTo("other"));
        lsensor = LocationClass.scheme(LocationClass.Type.STOP, false);
        lsensor.setLocs(newList);
        assertThat(lsensor.classOf(150), equalTo("other"));
        assertThat(lsensor.classOf(9999), equalTo("other"));
        assertThat(lsensor.classOf(8000), equalTo("other"));
        assertThat(lsensor.classOf(8597), equalTo("stop"));
        assertThat(lsensor.classOf(9502), equalTo("other"));
        lsensor = LocationClass.scheme(LocationClass.Type.STOP, true);
        lsensor.setLocs(newList);
        assertThat(lsensor.classOf(150), equalTo("other"));
        assertThat(lsensor.classOf(9999), equalTo("other"));
        assertThat(lsensor.classOf(8000), equalTo("other"));
        assertThat(lsensor.classOf(8597), equalTo("stop"));
        assertThat(lsensor.classOf(9502), equalTo("stop"));
    }

    /**
     * test codon filters
     */
    public void testCodonFilter() {
        String sequence = "AATGTGACCTGAATAATAG";
        CodonFilter filter = new CodonFilter("ATG", "AAT");
        assertTrue(filter.matches(1, sequence));
        assertTrue(filter.matches(2, sequence));
        assertFalse(filter.matches(3, sequence));
        assertFalse(filter.matches(4, sequence));
        assertFalse(filter.matches(5, sequence));
        assertFalse(filter.matches(6, sequence));
        assertFalse(filter.matches(7, sequence));
        assertTrue(filter.matches(12, sequence));
    }
}
