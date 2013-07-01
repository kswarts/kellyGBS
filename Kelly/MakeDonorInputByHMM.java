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
    public static byte diploidN= (byte) 0xff;
    private static double[][] trans;//transition matrix, based on arbitrary f in constructor
    private static double[][] emiss;//emission matrix, based on the probability of observing a state given the read depth and HW expected frequency for the site
    
    //*Takes in an hdf5 file, pulls out the heterozygous taxa, extract a matrix with depth info,*//
    //*calculate probabilities that actually homo/het for remaining sites based on HW expected *//
    //*frequencies and read depth. Does not take LD into account (assume independance for each site with respect to homo/het). *//
    //*Calculate joint probabilities based on 2pq for each site. Parse out the homozygous segments,*//
    //*returns the inbreds and inbred segments*//
    public static void newHomoSegDonorPrep(String inFileRoot, double hetCutoff, double f) {
        trans= new double[][] {
            {f,1-f},//homo {homo,het} //maybe add distance from
            {f,1-f}//het {homo,het}
        };
        MutableNucleotideAlignmentHDF5 mnah5= MutableNucleotideAlignmentHDF5.getInstance(dir+inFileRoot+".hmp.h5");
        System.out.println("read in h5 file");
        for (int taxon = 0; taxon < mnah5.getSequenceCount(); taxon++) {
            if (((double)mnah5.getHeterozygousCountForTaxon(taxon)/(double)mnah5.getTotalNotMissingForTaxon(taxon))<hetCutoff) continue;//if an inbred line, don't do anything
            System.out.println("Working on taxa"+mnah5.getTaxaName(taxon));
            double probTrueHomo= -1;
            double probTrueHet= -1;
            for (int site = 0; site < mnah5.getSequenceCount(); site++) {
                if (mnah5.getBase(taxon, site)==diploidN) continue; //skip sites with no data
                int obs= mnah5.isHeterozygous(taxon, site)==true?1:0;
                double currRD= getReadDepthForMinMaj(mnah5.getDepthForAlleles(taxon, site),
                        getIndexForMinMaj(mnah5.getMajorAllele(site),mnah5.getMinorAllele(site)));
                double currPQ= 2*mnah5.getMajorAlleleFrequency(site)*mnah5.getMinorAlleleFrequency(site);
                emiss= new double[][] {//initialize new emission matrix
                    {.999,.001},//homo {obsHomo,obsHet}
                    {currRD>1?((1-currPQ*currPQ*Math.pow(.5,currRD))/(currPQ*Math.pow(.5,currRD))):1,
                        currRD>1?((currPQ*currPQ*(1-Math.pow(.5,currRD)))/(currPQ*(1-Math.pow(.5,currRD)))):0}//het {obsHomo,obsHet}
                };
                if (probTrueHomo==-1) { //initialize HMM with HW expected allele freq
                    probTrueHomo= (1-currPQ)*emiss[0][obs]; //the starting probabilities for each true state, given the observed state
                    probTrueHet= currPQ*emiss[1][obs];
                    if (probTrueHet>probTrueHomo) mnah5.setBase(taxon, site, diploidN); //if true state chosen to be het, change site to missing
                    continue;
                }
                double homoHomo= probTrueHomo*trans[0][0]*emiss[0][obs]; //path probabilities from true homo at last site to true homo at curr site
                double homoHet= probTrueHomo*trans[0][1]*emiss[1][obs];
                double hetHomo= probTrueHet*trans[1][0]*emiss[0][obs];
                double hetHet= probTrueHet*trans[1][1]*emiss[1][obs];
                probTrueHomo=hetHomo>homoHomo?hetHomo:homoHomo;
                probTrueHet= hetHet>homoHet?hetHet:homoHet;
                if (probTrueHet>probTrueHomo) mnah5.setBase(taxon, site, diploidN);
            }
        }
        mnah5.clean();
        ExportUtils.writeToHDF5(mnah5, dir+inFileRoot+"inbredByHMM.hmp.h5");
    }
    
    private static void setToMissing(MutableNucleotideAlignmentHDF5 mnah5, int firstSite, int lastSite) {
        for (int site = firstSite; site < lastSite+1; site++) {
            mnah5.setBase(site, site, diploidN);
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
    
    public static void main(String[] args) {
        dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
        String fileName= "AllZeaGBSv27";
        newHomoSegDonorPrep(fileName,.05,.25);
    }
}
