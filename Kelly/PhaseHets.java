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
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableSingleEncodeAlignment;

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
    
    public static void phaseByHapFile(String hapFile, String fileToPhase, double mismatchRate) {
        Alignment haps= ImportUtils.readFromHapmap(hapFile, null);
        Alignment phase= ImportUtils.readFromHapmap(fileToPhase,null);
        MutableNucleotideAlignment newOne= MutableNucleotideAlignment.getInstance(phase.getIdGroup(), phase.getSiteCount());
        MutableNucleotideAlignment newTwo= MutableNucleotideAlignment.getInstance(phase.getIdGroup(), phase.getSiteCount());
        if (haps.getSiteCount()!=phase.getSiteCount()) System.out.println("Sites do not match between alignment and haplotypes");
        //based on 20 minor allele segments, choose top two haplotypes that best describe alignment segment. Go through and assign het alignment to new
        //alignment strands. Make it a sliding window so there is overlap and everything gets on the correct strand. If only one looks good for a certain
        //window, just go based on that
        for (int site=0;site<phase.getSequenceCount();site++) {
            
        }
    }
    
    public static void getClosestHaps() {
        
    }
    
    public static void GetExpectedHomozygosity(String inFile, double minCov, double minMAF, double hetCutoff, boolean gz) {
        Alignment a= ImportUtils.readFromHapmap(dir+inFile+(gz==true?".hmp.txt.gz":".hmp.txt"), null);
        System.out.println("reading from file: "+dir+inFile+(gz==true?".hmp.txt.gz":".hmp.txt"));
        //filter for MAF
        a= KellyUtils.SubsetHapmapByMAF(a, minMAF, false);
        //filter for coverage\
        a= KellyUtils.SubsetHapmapByTaxaCov(a, minCov, false, true);
        Alignment highHet= KellyUtils.SubsetHapmapByHeterozygosity(a, hetCutoff, true, false,true);
        Alignment lowHet= KellyUtils.SubsetHapmapByHeterozygosity(a, hetCutoff, false, false,true);
        if (highHet.getSequenceCount()>0) {
            ExportUtils.writeToHapmap(highHet, true, dir+inFile+"highHet.hmp.txt.gz", '\t', null);
            Alignment newHighHet= ImportUtils.readFromHapmap(dir+inFile+"highHet.hmp.txt.gz", null);
            ArrayList<Double> highHetHomoSeg= new ArrayList<Double>(newHighHet.getSiteCount()*newHighHet.getSequenceCount());
            for (int taxon=0; taxon<newHighHet.getSequenceCount(); taxon++) {
                int highCount= 0;
                for (int site= 0; site<newHighHet.getSiteCount(); taxon++) {
                    if (newHighHet.getBaseArray(taxon, site)[0]==newHighHet.getBaseArray(taxon, site)[1]) highCount++;
                    else {
                        highHetHomoSeg.add(Math.pow(1-minMAF,highCount));
                        highCount= 0;
                    }
                }
            }
            highHetHomoSeg.trimToSize();
            
            try{
            DataOutputStream highHetsOut= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(dir+inFile+"_highHets.txt"), 655360));
            highHetsOut.writeBytes("in file name: "+inFile+"\nminCov: "+minCov+"\nminMAF: "+minMAF+"\nHet cutoff: "+hetCutoff);
            for (int high= 0; high < highHetHomoSeg.size(); high++) {
                highHetsOut.writeBytes("\n"+highHetHomoSeg.get(high));
            }
            highHetsOut.close();
            }

           catch(IOException e) {
                System.out.println(e);
            }
        }
        if (lowHet.getSequenceCount()>0) {
            ExportUtils.writeToHapmap(lowHet, true, dir+inFile+"lowHet.hmp.txt.gz", '\t', null);
            Alignment newLowHet= ImportUtils.readFromHapmap(dir+inFile+"lowHet.hmp.txt.gz", null);
            ArrayList<Double> lowHetHomoSeg= new ArrayList<Double>(newLowHet.getSiteCount()*newLowHet.getSequenceCount());
            for (int taxon=0; taxon<newLowHet.getSequenceCount(); taxon++) {
                int lowCount= 0;
                for (int site= 0; site<newLowHet.getSiteCount(); taxon++) {
                    if (newLowHet.getBaseArray(taxon, site)[0]==newLowHet.getBaseArray(taxon, site)[1]) lowCount++;
                    else {
                        lowHetHomoSeg.add(Math.pow(1-minMAF,lowCount));
                        lowCount= 0;
                    }
                }
            }
            lowHetHomoSeg.trimToSize();
            
            try{
            DataOutputStream lowHetsOut= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(dir+inFile+"_lowHets.txt"), 655360));
            lowHetsOut.writeBytes("in file name: "+inFile+"\nminCov: "+minCov+"\nminMAF: "+minMAF+"\nHet cutoff: "+hetCutoff);
            for (int low= 0; low < lowHetHomoSeg.size(); low++) {
                lowHetsOut.writeBytes("\n"+lowHetHomoSeg.get(low));
            }
            lowHetsOut.close();
            }

           catch(IOException e) {
                System.out.println(e);
            }
        }
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
        
        GetExpectedHomozygosity("AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_minHet0.024",.7,.1,.15,true);
    }
}
