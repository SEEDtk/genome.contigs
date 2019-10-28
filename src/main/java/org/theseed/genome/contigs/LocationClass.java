/**
 *
 */
package org.theseed.genome.contigs;

import org.theseed.locations.Frame;
import org.theseed.locations.LocationList;

/**
 * This is the base class for location classifiers.  It computes the class of a location
 * (usually related to its coding function) based on the current and previous frame
 * codes.  The constructor specifies whether or not minus strand proteins are considered
 * as coding regions.
 *
 * @author Bruce Parrello
 *
 */
public abstract class LocationClass {


    // FIELDS
    /** If TRUE, then minus-strand proteins are considered coding regions */
    boolean	negative;
    /** controlling location list */
    LocationList contigLocs;


    /**
     * Construct a blank location class handler.
     *
     * @param negativeFlag	TRUE if minus-strand proteins are considered coding regions
     */
    public LocationClass(boolean negativeFlag) {
        this.negative = negativeFlag;
        this.contigLocs = null;
    }

    /**
     * Return the location class.  A return of NULL means the current location
     * is invalid.
     *
     * @param pos	position of the location whose class is desired
     */
    public abstract String classOf(int pos);

    /**
     * Convert a frame according to the policy on the minus strand.  This
     * will either return the original frame or Frame.F0.
     *
     * @param frm	frame in question
     */
    public Frame normalize(Frame frm) {
        Frame retVal = frm;
        if (! this.negative && frm.negative()) {
            retVal = Frame.F0;
        }
        return retVal;
    }

    /** enum used to choose a location class scheme */
    public static enum Type {
        /** identify "coding" or "space" (non-coding) locations */
        CODING,
        /** identify "start" (start of coding), "stop" (past end of coding), or "other" locations */
        EDGE,
        /** identify actual coding frame */
        PHASE,
        /** identify "start" (start of coding) or "other" locations */
        START,
        /** identify "stop" (past end of coding) or "other" locations */
        STOP
    }

    /**
     * @return the appropriate codon filter for the location class, or NULL if none is needed
     *
     * @param type	location class scheme
     */
    public static CodonFilter filter(Type type) {
        CodonFilter retVal;
        switch (type) {
        case EDGE :
            retVal = new CodonFilter("ATG", "GTG", "TTG", "TAA", "TAG", "TGA");
            break;
        case START :
            retVal = new CodonFilter("ATG", "GTG", "TTG");
            break;
        case STOP :
            retVal = new CodonFilter("TAA", "TAG", "TGA");
            break;
        default :
            retVal = null;
        }
        return retVal;
    }

    public static LocationClass scheme(Type type, boolean negativeFlag) {
        LocationClass retVal;
        switch (type) {
        case CODING :
            retVal = new LocationClass.Coding(negativeFlag);
            break;
        case EDGE :
            retVal = new LocationClass.Edge(negativeFlag);
            break;
        case PHASE :
            retVal = new LocationClass.Phase(negativeFlag);
            break;
        case START :
            retVal = new LocationClass.Start(negativeFlag);
            break;
        case STOP :
            retVal = new LocationClass.Stop(negativeFlag);
            break;
        default :
            retVal = null;
        }
        return retVal;
    }

    /** Store the controlling location list */
    public void setLocs(LocationList contigLocs) {
        this.contigLocs = contigLocs;
    }

    // SUBCLASSES

    /**
     * Classification is the actual coding frame string (+1, +2, +3, 0)
     */
    public static class Phase extends LocationClass {

        public Phase(boolean negativeFlag) {
            super(negativeFlag);
        }

        @Override
        public String classOf(int pos) {
            String retVal = null;
            Frame frm = this.contigLocs.computeRegionFrame(pos, pos);
            if (frm != Frame.XX)
                retVal = this.normalize(frm).toString();
            return retVal;
        }

    }

    /**
     * Classification is "start", "stop", or "other".
     */
    public static class Edge extends LocationClass {

        public Edge(boolean negativeFlag) {
            super(negativeFlag);
        }

        @Override
        public String classOf(int pos) {
            LocationList.Edge type = this.contigLocs.isEdge(pos, this.negative);
            return type.toString();
        }


    }

    /**
     * Classification is "coding" or "space".
     */
    public static class Coding extends LocationClass {

        public Coding(boolean negativeFlag) {
            super(negativeFlag);
        }

        @Override
        public String classOf(int pos) {
            Frame frm = this.contigLocs.computeRegionFrame(pos, pos);
            frm = this.normalize(frm);
            String retVal;
            switch (frm) {
            case XX :
                retVal = null;
                break;
            case F0 :
                retVal = "space";
                break;
            default:
                retVal = "coding";
            }
            return retVal;
        }

    }

    /**
     * Classification is "stop" or "other".
     */
    public static class Stop extends LocationClass {

        public Stop(boolean negativeFlag) {
            super(negativeFlag);
        }

        @Override
        public String classOf(int pos) {
            LocationList.Edge type = this.contigLocs.isEdge(pos, this.negative);
            return (type == LocationList.Edge.STOP ? "stop" : "other");
        }

    }

    /**
     * Classification is "start" or "other".
     */
    public static class Start extends LocationClass {

        public Start(boolean negativeFlag) {
            super(negativeFlag);
        }

        @Override
        public String classOf(int pos) {
            LocationList.Edge type = this.contigLocs.isEdge(pos, this.negative);
            return (type == LocationList.Edge.START ? "start" : "other");
        }

    }

}
