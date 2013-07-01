/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 *
 * @author kelly
 */
public class MakeDonorInputByHMM {
    public static String dir;
    private static boolean[] hetTaxa;
    private static int[][] covMask;//[taxon][site]
    private static double currRD;
    private static double[] twoPQ;//the HW expectations for each site
    private static double currPQ;
    public static byte diploidN= (byte) 0xff;
    private static double inbreedCoef;
    private static double[] start= new double[] {.2,.8};//starting probabilities for {homo/het}
    private static double[][] trans= new double[][] {
        {inbreedCoef,1-inbreedCoef},//homo {homo,het} //maybe add distance from
        {inbreedCoef,1-inbreedCoef}//het {homo,het}
    };
    //based on the probability of observing a state given the read depth and HW expected frequency for the site
    private static double[][] emiss= new double[][] {
        {.999,.001},//homo {obsHomo,obsHet}
        {currRD>1?((1-currPQ*currPQ*Math.pow(.5,currRD))/(currPQ*Math.pow(.5,currRD))):1,
            currRD>1?((currPQ*currPQ*(1-Math.pow(.5,currRD)))/(currPQ*(1-Math.pow(.5,currRD)))):0}//het {obsHomo,obsHet}
    };
    
    //*Takes in an hdf5 file, pulls out the heterozygous taxa, extract a matrix with depth info,*//
    //*calculate probabilities that actually homo/het for remaining sites based on HW expected *//
    //*frequencies and read depth. Does not take LD into account (assume independance for each site with respect to homo/het). *//
    //*Calculate joint probabilities based on 2pq for each site. Parse out the homozygous segments,*//
    //*returns the inbreds and inbred segments*//
    public static void newHomoSegDonorPrep(String inFileRoot, double hetCutoff, double f) {
        inbreedCoef= f;
        MutableNucleotideAlignmentHDF5 mnah5= MutableNucleotideAlignmentHDF5.getInstance(dir+inFileRoot+".imp.hmp.h5");
        getHetTaxa(mnah5,hetCutoff);
        getMaskForReadDepth(mnah5);
        setProbHet(mnah5);
        for (int taxon = 0; taxon < mnah5.getSequenceCount(); taxon++) {
            if (hetTaxa[taxon]==false) continue;//if an inbred line, don't do anything
            boolean started= false;
            double hetHet;//the path probability from het to het
            double hetHomo;//the path probability from het to homo
            double homoHomo;
            double homoHet;
            double probGivenHomo= 0;
            double probGivenHet= 0;
            for (int site = 0; site < mnah5.getSequenceCount(); site++) {
                if (mnah5.getBase(taxon, site)==diploidN) continue; //skip sites with no data
                int obs= mnah5.isHeterozygous(taxon, site)==true?1:0;
                currRD= covMask[taxon][site];
                currPQ= twoPQ[site];
                if (started==false) {
                    probGivenHomo= start[0]*emiss[0][obs]; //the starting probabilities for each true state, given the observed state
                    probGivenHet= start[1]*emiss[1][obs];
                    if (probGivenHet>probGivenHomo) mnah5.setBase(taxon, site, diploidN); //if true state chosen to be het, change site to missing
                    started= true;
                    continue;
                }
                homoHomo= probGivenHomo*trans[0][0]*emiss[0][obs];
                homoHet= probGivenHomo*trans[0][1]*emiss[1][obs];
                hetHomo= probGivenHet*trans[1][0]*emiss[0][obs];
                hetHet= probGivenHet*trans[1][1]*emiss[1][obs];
                
                probGivenHomo=hetHomo>homoHomo?hetHomo:homoHomo;
                probGivenHet= hetHet>hetHomo?hetHet:hetHomo;
                if (probGivenHet>probGivenHomo) mnah5.setBase(taxon, site, diploidN);
            }
        }
        mnah5.clean();
        ExportUtils.writeToHDF5(mnah5, dir+inFileRoot+"inbredByHW.imp.hmp.h5");
    }
    
    private static void setToMissing(MutableNucleotideAlignmentHDF5 mnah5, int firstSite, int lastSite) {
        for (int site = firstSite; site < lastSite+1; site++) {
            mnah5.setBase(site, site, diploidN);
        }
    }
    
    private static void getHetTaxa(MutableNucleotideAlignmentHDF5 mnah5, double hetCutoff) {
        hetTaxa= new boolean[mnah5.getSequenceCount()];
        for (int taxon = 0; taxon < mnah5.getSequenceCount(); taxon++) {
            if (((double)mnah5.getHeterozygousCountForTaxon(taxon)/(double)mnah5.getTotalNotMissingForTaxon(taxon))>hetCutoff) hetTaxa[taxon]= true;
        }
    }
    
    private static void getMaskForReadDepth(MutableNucleotideAlignmentHDF5 mnah5) {
        covMask= new int[mnah5.getSequenceCount()][mnah5.getSiteCount()];
        for (int taxon = 0; taxon < mnah5.getSequenceCount(); taxon++) {
            for (int site = 0; site < mnah5.getSiteCount(); site++) {
                covMask[taxon][site]= getReadDepthForMinMaj(mnah5.getDepthForAlleles(taxon, site),
                        getIndexForMinMaj(mnah5.getMajorAllele(site),mnah5.getMinorAllele(site)));
            }
        }
    }
    
    private static int[] getIndexForMinMaj(byte maj, byte min) {
        int[] minMaj= new int[2];
        minMaj[0]= maj==NucleotideAlignmentConstants.A_ALLELE?0:maj==NucleotideAlignmentConstants.C_ALLELE?
                1:maj==NucleotideAlignmentConstants.G_ALLELE?2:maj==NucleotideAlignmentConstants.T_ALLELE?
                3:maj==NucleotideAlignmentConstants.INSERT_ALLELE?4:5;
        minMaj[1]= min==NucleotideAlignmentConstants.A_ALLELE?0:min==NucleotideAlignmentConstants.C_ALLELE?
                1:min==NucleotideAlignmentConstants.G_ALLELE?2:min==NucleotideAlignmentConstants.T_ALLELE?
                3:min==NucleotideAlignmentConstants.INSERT_ALLELE?4:5;
        return minMaj;
    }
    
    public static int getTotalReadDepth(byte[] depthForAlleles) {
        int depth= 0;
        for (int all = 0; all < depthForAlleles.length; all++) {
            depth+= depthForAlleles[all];
        }
        return depth;
    }
    
    public static int getReadDepthForMinMaj(byte[] depthForAlleles, int[]minMajIndices) {
        int depth= (int)depthForAlleles[minMajIndices[0]]+(int)depthForAlleles[minMajIndices[1]];
        return depth;
    }
    
    private static void setProbHet(MutableNucleotideAlignmentHDF5 mnah5) {
        twoPQ= new double[mnah5.getSiteCount()];
        for (int site = 0; site < mnah5.getSiteCount(); site++) {
            twoPQ[site]= 2*mnah5.getMajorAlleleFrequency(site)*mnah5.getMinorAlleleFrequency(site);
        }
    }
    
    public static void main(String[] args) {
        dir= "/Users/kelly/Documents/GBS/Imputation/SmallFiles/";
        String fileName= "all26.130530.ALL";
        newHomoSegDonorPrep(fileName,.15,.25);
    }
}
