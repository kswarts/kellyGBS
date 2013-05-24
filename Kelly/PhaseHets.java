/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;

/**
 *
 * @author kls283
 */
public class PhaseHets {
    //Run Ed's code to get haplotypes for whole dataset
    //Only keep those with no hets
    //Use haplotypes to phase taxa (store two strands in different alignments)
    //Markov chain given nearest site as prior state. start with first present allele on one alignment, add alleles to alignment or alternate alignment based on conditional probability of seeing the states together,based on haplotypes
    public static final byte N= (byte) 0xFF;
    public static final byte A = (byte) 0x00;
    public static final byte C = (byte) 0x11;
    public static final byte G = (byte) 0x22;
    public static final byte T = (byte) 0x33;
    public static final byte INSERT = (byte) 0x44;
    public static final byte GAP = (byte) 0x55;
    
    public static String dir;
    public static double[] dist;
    
    public static void phaseByHapFile(String hapFile, String fileToPhase, double mismatchTol, int windowSize) {
        Alignment haps= ImportUtils.readFromHapmap(hapFile, null);
        Alignment phase= ImportUtils.readFromHapmap(fileToPhase,null);
        MutableNucleotideAlignment newOne= MutableNucleotideAlignment.getInstance(phase.getIdGroup(), phase.getSiteCount());
        MutableNucleotideAlignment newTwo= MutableNucleotideAlignment.getInstance(phase.getIdGroup(), phase.getSiteCount());
        if (haps.getSiteCount()!=phase.getSiteCount()) System.out.println("Sites do not match between alignment and haplotypes");
        //based on 20 minor allele segments, choose top two haplotypes that best describe alignment segment. Go through and assign het alignment to new
        //alignment strands. Make it a sliding window so there is overlap and everything gets on the correct strand. If only one looks good for a certain
        //window, just go based on that
        for (int window= 0;window<(phase.getSiteCount()/windowSize)+1;window++) {
            FilterAlignment hapsForWindow= FilterAlignment.getInstance(haps, window, (window+windowSize<phase.getSiteCount())?window+windowSize:phase.getSiteCount());
            FilterAlignment taxonForWindow= FilterAlignment.getInstance(phase, window, (window+windowSize<phase.getSiteCount())?window+windowSize:phase.getSiteCount());
            
            byte[] diploidMinor= new byte[taxonForWindow.getSiteCount()]; //make an array of diploid minor allele states for sites in window
            for (int site = 0; site < taxonForWindow.getSiteCount(); site++) {
                    diploidMinor[site]= getDiploidBase(taxonForWindow.getMinorAllele(site));
                }
            
            for (int taxon=0;taxon<phase.getSequenceCount();taxon++) {
                int[] bestHaps= getExplanatoryHaps(hapsForWindow,taxonForWindow,taxon);//get the 10 most similar haplotypes by distance, ranked by index in haps
                
                //check to make sure that top haplotypes explain the minor alleles, up to tolerated mismatch
                int secondTaxon= -1;
                for (int hapMatch = 1; hapMatch < bestHaps.length; hapMatch++) {
                    if (secondTaxon>0) break;
                    int mismatches= 0;
                    for (int site= 0;site<phase.getSiteCount();site++) {
                        if (taxonForWindow.getBase(taxon, site)!=diploidMinor[site]) continue;
                        if(taxonForWindow.getBase(taxon, site)==hapsForWindow.getBase(bestHaps[0], site)||
                                taxonForWindow.getBase(taxon, site)==hapsForWindow.getBase(bestHaps[hapMatch], site)) continue;
                        else mismatches++;
                        if (mismatches>mismatchTol) break;
                        secondTaxon= hapMatch;
                    }
                }
                //use the top two haplotypes to phase alignment   
            }
        }
    }
    
    public static int[] getExplanatoryHaps(FilterAlignment localHaps, FilterAlignment localPhase, int taxonToPhase) {
        long[] pMj= localPhase.getAllelePresenceForAllSites(taxonToPhase, 0).getBits();
        long[] pMn= localPhase.getAllelePresenceForAllSites(taxonToPhase, 1).getBits();
        for (int hap=0;hap<localHaps.getSequenceCount();hap++) {
            long[] hMj= localHaps.getAllelePresenceForAllSites(taxonToPhase, 0).getBits();
            long[] hMn= localHaps.getAllelePresenceForAllSites(taxonToPhase, 1).getBits();
            dist[hap]= IBSDistanceMatrix.computeHetBitDistances(hMj, hMn, pMj, pMn, 50)[0];
        }
        double[] sortDist= dist;
        Arrays.sort(sortDist, 0, dist.length-1);
        int[] bestHaps= new int[10];
        for (int index = 0; index < 10; index++) {
            bestHaps[index]= Arrays.binarySearch(dist, sortDist[index]);
            
        }
        return bestHaps;
    }
    
    public static void phaseByMarkov(String hapFile, String fileToPhase, String outFile) {
        Alignment haps= ImportUtils.readFromHapmap(hapFile, null);
        Alignment phase= ImportUtils.readFromHapmap(fileToPhase,null);
        if (haps.getSiteCount()!=phase.getSiteCount()) System.out.println("Sites do not match between alignment and haplotypes");
        Alignment[] doubleAlign= {phase,phase};
        MutableSingleEncodeAlignment gametes= MutableSingleEncodeAlignment.getInstance(doubleAlign);
        int plus= phase.getSequenceCount();
        int lastSite= 0;
        for (int taxon= 0; taxon<phase.getSequenceCount();taxon++) {
            for (int site= 0; site<phase.getSiteCount(); site++) {
                if (phase.getBase(taxon, site)!=N) {
                    if (lastSite==0) {//for the first not missing site the algorithm comes across
                        if (phase.isHeterozygous(taxon, site)==true) {//if het fill in appropriate diploid bases for both versions of the alignment
                            gametes.setBase(taxon, site, PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[0]));
                            gametes.setBase(taxon+plus, site, PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[1]));
                        }
                        else gametes.setBase(taxon+plus, site, N);//set second alignment to missing
                        lastSite= site;
                    }
                    else {
                        if (gametes.getBase(taxon+plus, lastSite)==N) {//if only working on first alignment
                            if (phase.isHeterozygous(taxon, site)==true) {
                                byte firstBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[0]);
                                byte otherBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[1]);
                                if (Math.random() <= getConditional(haps,site,firstBase,lastSite,phase.getBase(taxon, lastSite))) {
                                    gametes.setBase(taxon, site, firstBase);
                                    gametes.setBase(taxon+plus, site, otherBase);
                                }
                                else {
                                    gametes.setBase(taxon, site, otherBase);
                                    gametes.setBase(taxon+plus, site, firstBase);
                                }
                            }
                            else {
                                if (Math.random() <= getConditional(haps,site,phase.getBase(taxon, site),lastSite,phase.getBase(taxon, lastSite))) gametes.setBase(taxon+plus, site, N);
                                else gametes.setBase(taxon, site, N);
                            }
                        }
                        else if (gametes.getBase(taxon, lastSite)==N) {//if only working on second alignment
                            if (phase.isHeterozygous(taxon, site)==true) {
                                byte firstBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[0]);
                                byte otherBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[1]);
                                if (Math.random() <= getConditional(haps,site,firstBase,lastSite,phase.getBase(taxon, lastSite))) {
                                    gametes.setBase(taxon+plus, site, firstBase);
                                    gametes.setBase(taxon, site, otherBase);
                                }
                                else {
                                    gametes.setBase(taxon+plus, site, otherBase);
                                    gametes.setBase(taxon, site, firstBase);
                                }
                            }
                            else {
                                if (Math.random() <= getConditional(haps,site,phase.getBase(taxon, site),lastSite,phase.getBase(taxon, lastSite))) gametes.setBase(taxon, site, N);
                                else gametes.setBase(taxon+plus, site, N);
                            }
                        }
                        else {//if both alignments have alleles present at the previous site
                            if (phase.isHeterozygous(taxon, site)==true) {
                                byte firstBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[0]);
                                byte otherBase= PhaseHets.getDiploidBase(phase.getBaseArray(taxon, site)[1]);
                                double probFirstAlignOne= getConditional(haps,site,firstBase,lastSite,gametes.getBase(taxon, lastSite));
                                double probFirstAlignTwo= getConditional(haps,site,firstBase,lastSite,gametes.getBase(taxon+plus, lastSite));
                                if (probFirstAlignOne>=probFirstAlignTwo) {//assign the first base to the alignment with the highest cond probability
                                    gametes.setBase(taxon, site, firstBase);
                                    gametes.setBase(taxon+plus, site, otherBase);
                                }
                                else {
                                    gametes.setBase(taxon+plus, site, firstBase);
                                    gametes.setBase(taxon, site, otherBase);
                                }
                            }
                            else {
                                byte currBase= phase.getBase(taxon, site);
                                double probAlignOne= getConditional(haps,site,currBase,lastSite,gametes.getBase(taxon, lastSite));
                                double probAlignTwo= getConditional(haps,site,currBase,lastSite,gametes.getBase(taxon+plus, lastSite));
                                if (probAlignOne>=probAlignTwo) {//assign the base to the alignment with the highest conditional probability and set the other alignment to missing at this site
                                    gametes.setBase(taxon, site, currBase);
                                    gametes.setBase(taxon+plus, site, N);
                                }
                                else {
                                    gametes.setBase(taxon+plus, site, currBase);
                                    gametes.setBase(taxon, site, N);
                                }
                            }
                        }
                        lastSite= site;
                    }
                }
            }
        }
        gametes.clean();
        ExportUtils.writeToHapmap(gametes, true, outFile, '\t', null);
    }
    
    public static double getConditional(Alignment haps, int focalSite, byte focalBase, int conditionalSite, byte conditionalBase) {
        double numCond= 0.0;
        double match= 0.0;
        for (int hap= 0;hap<haps.getSequenceCount();hap++) {
            if (haps.getBase(hap, conditionalSite)==conditionalBase) {
                numCond++;
                if (haps.getBase(hap, focalSite)==focalBase) match++;
            }
        }
        double cond= match/numCond!=0?match/numCond:0;
        
        return cond;
    }
    
    public static byte getDiploidBase(byte haploidBase) {
        byte diploidBase;
        if (haploidBase==(byte) 0xF) diploidBase= N;
        else if (haploidBase==(byte) 0x0) diploidBase= A;
        else if (haploidBase==(byte) 0x1) diploidBase= C;
        else if (haploidBase==(byte) 0x2) diploidBase= G;
        else if (haploidBase==(byte) 0x3) diploidBase= T;
        else if (haploidBase==(byte) 0x4) diploidBase= GAP;
        else diploidBase= INSERT;
        return diploidBase;
    }
    
    public static void main(String[] args) {
        dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
        String hapmap= dir+"AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1.hmp.txt";
//        String haplotype= dir+"AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_HaplotypeMerge.hmp.txt";
//        String outFile= dir+"SEED_12S_GBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_PhasedFromAllZeaHaplotypeMerge.hmp.txt";
//        phaseByMarkov(hapmap,haplotype,outFile);
        
    }
}
