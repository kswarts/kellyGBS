/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.FileInputStream;
import java.util.*;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;

/**
 *
 * @author kelly
 */
public class MakeDonorInputByHMM {
    public static String dir;
    public static byte diploidN= NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
    private static double[][] trans;//transition matrix, based on arbitrary f in constructor
    private static double[][] emiss;//emission matrix, based on the probability of observing a state given the read depth and HW expected frequency for the site
    private static MutableNucleotideAlignmentHDF5 mnah5;
    private static int[] whichWindow; //length of getSiteCount, hold window
    private static HashMap<Integer,Double> rho; //window->rho
    private static int[] windowOffset; //window->chromosome
    private static String[] whichChr;
    //*Takes in an hdf5 file, pulls out the heterozygous taxa, extract a matrix with depth info,*//
    //*calculate probabilities that actually homo/het for remaining sites based on HW expected *//
    //*frequencies and read depth. Does not take LD into account (assume independance for each site with respect to homo/het). *//
    //*Calculate joint probabilities based on 2pq for each site. Parse out the homozygous segments,*//
    //*returns the inbreds and inbred segments*//
    
    public static void newHomoSegDonorPrep(String inFile, String recombFile, double hetCutoff, double f) {//this one only acts on LD and readDepth
        mnah5= MutableNucleotideAlignmentHDF5.getInstance(dir+inFile);
        System.out.println("read in h5 file");
        readRecombinationFromTabTxt(recombFile, 10);
        Alignment a= ImportUtils.readGuessFormat(dir+inFile, false);
        MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
        System.out.println("generate mna");
        if (mnah5.getLocus(0).getChromosomeName().equalsIgnoreCase(whichChr[0])==false||mnah5.getLocus(mnah5.getSiteCount()-1).getChromosomeName().equalsIgnoreCase(whichChr[whichChr.length-1])==false) 
            System.out.println("Chromosomes do not match between hapmap and recombination file. Sites with missing rho will be skipped");
        for (int taxon = 0; taxon < mnah5.getSequenceCount(); taxon++) {
            if (((double)mnah5.getHeterozygousCountForTaxon(taxon)/(double)mnah5.getTotalNotMissingForTaxon(taxon))<hetCutoff) continue;//if an inbred line, don't do anything
            System.out.println("Working on taxon"+mnah5.getTaxaName(taxon));
            double probTrueHomo= -1;
            double probTrueHet= -1;
            boolean newChr= true;
            int lastSite= 0;
            for (int site = 0; site < mnah5.getSequenceCount(); site++) {
                if (mnah5.getLocusName(site).equalsIgnoreCase(whichChr[whichChr.length-1])) continue;//skip chromosomes with no recombination information
                if (mnah5.getMinorAlleleFrequency(site)==0) continue; //skip invariant sites
                if (mnah5.getBase(taxon, site)==diploidN) continue; //skip sites with no data
                if (mnah5.getLocus(site).getChromosomeName().equalsIgnoreCase(mnah5.getLocus(lastSite).getChromosomeName())==false) newChr= true;
                int obs= mnah5.isHeterozygous(taxon, site)==true?1:0;
                double currRD= getReadDepthForMinMaj(mnah5.getDepthForAlleles(taxon, site),
                        getIndexForMinMaj(mnah5.getMajorAllele(site),mnah5.getMinorAllele(site)));
                double currPQ= 2*mnah5.getMajorAlleleFrequency(site)*mnah5.getMinorAlleleFrequency(site);
//                emiss= new double[][] {//initialize new emission matrix using just read depth
//                    {.999,.001},//homo {obsHomo,obsHet}
//                    {currRD>1?Math.pow(.5,currRD):.999,
//                        currRD>1?1.0-Math.pow(.5,currRD):0.001}//het {obsHomo,obsHet}
//                };
                emiss= new double[][] {//initialize new emission matrix using HW/RD joint probability
                    {.999,.001},//homo {obsHomo,obsHet}
                    {currRD>1?(((1.0-currPQ)*currPQ*Math.pow(.5,currRD))/(currPQ*Math.pow(.5,currRD))):.99999,
                        currRD>1?((currPQ*currPQ*(1.0-Math.pow(.5,currRD)))/(currPQ*(1.0-Math.pow(.5,currRD)))):0.00001}//het {obsHomo,obsHet}
                };
                if (newChr==true) { //initialize HMM with HW expected allele freq
                    probTrueHomo= (1-currPQ)*emiss[0][obs]; //the starting probabilities for each true state, given the observed state
                    probTrueHet= currPQ*emiss[1][obs];
                    if (probTrueHet>probTrueHomo) mnah5.setBase(taxon, site, diploidN); //if true state chosen to be het, change site to missing
                    lastSite= site; //this will result in strange LD for the next site if this is a het, but it only happens once and viterbi is bad here anyway
                    newChr= false;
                    continue;
                }
                
                if (mnah5.isHeterozygous(taxon, site)==true) {
                    probTrueHomo=0.001;
                    probTrueHet= .999;
                    lastSite= site;
                }
                else {
//                    double r= getLD(taxon,lastSite,site);
//                    trans= new double[][] {
//                        {r,(1.0-r)},//homo {homo,het}
//                        {(1.0-r),r}//het {homo,het}
//                    };
                    //the joint probability of recombination and inbreeding (f is made up)
                    double noRecomb= Math.pow((1-rho.get(whichWindow[site])),(mnah5.getPositionInLocus(site)-mnah5.getPositionInLocus(lastSite)));//prob of no recombination based on rho and the physical distance between the curr and last site
                    trans= new double[][] {
                        {noRecomb*f,(1-noRecomb)*(1-f)},//homo {homo,het}
                        {(1-noRecomb)*f,noRecomb*(1-f)}//het {homo,het}
                    };
                    double homoHomo= probTrueHomo*trans[0][0]*emiss[0][obs]; //path probabilities from true homo at last site to true homo at curr site
                    double homoHet= probTrueHomo*trans[0][1]*emiss[1][obs];
                    double hetHomo= probTrueHet*trans[1][0]*emiss[0][obs];
                    double hetHet= probTrueHet*trans[1][1]*emiss[1][obs];
                    probTrueHomo=hetHomo>homoHomo?hetHomo:homoHomo;
                    probTrueHet= hetHet>homoHet?hetHet:homoHet;
                    lastSite= site;
                }
                //if likely in a het sequence, set to missing
                if (probTrueHet>probTrueHomo) mna.setBase(taxon, site, diploidN);
            }
        }
        mna.clean();
        ExportUtils.writeToHapmap(mna, true, dir+inFile+"inbredByHMMWithHW.hmp.txt.gz", '\t', null);
    }
    
    private static double getLD(int taxon, int siteOne, int siteTwo) { //only use when target taxon is homozygous at both sites, assume both polymorphic
        byte siteOneAll= mnah5.getBaseArray(taxon, siteOne)[0];
        byte siteTwoAll= mnah5.getBaseArray(taxon, siteTwo)[0];
        double pOne= 0;
        double qOne= 0;
        double sitesTested= 0;
        double sitesMatch= 0;
        for (int t = 0; t < mnah5.getSequenceCount(); t++) { //get frequency of target haplotype
            boolean oneTrue= false;
            if (mnah5.getBase(taxon, siteOne)==diploidN||mnah5.getBase(taxon, siteTwo)==diploidN) continue;
            sitesTested++;
            if (mnah5.getBaseArray(t, siteOne)[0]==siteOneAll||mnah5.getBaseArray(t, siteOne)[1]==siteOneAll) {
                pOne++;
                oneTrue= true;
            }
            if (mnah5.getBaseArray(t, siteTwo)[0]==siteTwoAll||mnah5.getBaseArray(t, siteTwo)[1]==siteTwoAll) {
                qOne++;
                if (oneTrue==true) sitesMatch++;
            }
        }
        double siteOneFreq= pOne/sitesTested;
        double siteTwoFreq= qOne/sitesTested;
        double d= (sitesMatch/sitesTested)-(siteOneFreq*siteTwoFreq);
        double r= d<0?d/Math.min(siteOneFreq*siteTwoFreq, (1.0-siteOneFreq)*(1.0-siteTwoFreq)):d/Math.min(siteOneFreq*(1.0-siteTwoFreq), (1.0-siteOneFreq)*siteTwoFreq);
        return r;
    }
    
    //file should be tab delimited and organized as chromosome, starting position, estimated rho
    private static void readRecombinationFromTabTxt(String inFile, int numChr) {
        whichWindow= new int[mnah5.getSiteCount()];//hold the window based on site indices
        rho= new HashMap<Integer,Double>();//keys are the window, values are recomb rates
        windowOffset= new int[numChr];//the first window for each chromosome
        whichChr= new String[numChr];//an array of the chromosomes present in the recomb file
        ArrayList<String>[] physicalPos= new ArrayList[numChr];
        for (int i= 0;i<numChr;i++) {physicalPos[i]= new ArrayList<String>();}
        int index= 0;
        try {
            FileInputStream fis= new FileInputStream(inFile);
            Scanner scanner= new Scanner(fis);
            int currChr= 0;
            do {
                String[] next= scanner.nextLine().split("\t");
                int test= Integer.parseInt(next[0]);
                if (currChr!=Integer.parseInt(next[0])) {
                    windowOffset[currChr]= index;
                    whichChr[currChr]= next[0];
                } 
                currChr= Integer.parseInt(next[0]);
                physicalPos[currChr-1].add(next[1]);
                rho.put(index,Double.parseDouble(next[2]));
                index++;
                if (index==12593) {
                    int i= 0;
                }
            }
            while (scanner.hasNextLine());
            scanner.close();
            fis.close();
        }
        catch (Exception e) {
            System.out.println("Problem reading recombination file\nFile terminated at index "+index);
        }
        int[] startPosArray;
        for (int currChr = 0; currChr < physicalPos.length; currChr++) {
            startPosArray= new int[physicalPos[currChr].size()];
            for (int i= 0;i<physicalPos[currChr].size();i++) {startPosArray[i]= Integer.parseInt(physicalPos[currChr].get(i));}
            for (int site = mnah5.getLocus(whichChr[currChr]).getStart(); site < mnah5.getLocus(whichChr[currChr]).getEnd()+1; site++) {
                int currPos= mnah5.getPositionInLocus(site);
                int findWindow= Arrays.binarySearch(startPosArray, currPos);
                if (findWindow>-1) whichWindow[site]= windowOffset[currChr]+findWindow;
                else whichWindow[site]= windowOffset[currChr]+(-findWindow)-2;
            }
        }
        System.out.println("Recombination file read in with "+whichChr.length+" chromosomes and "+rho.size()+" windows");
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
//        dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
        dir= "/Users/kelly/Documents/GBS/Imputation/";
        String fileName= "AllZeaGBS_v2.7_SeqToGenos_part14.hmp.h5";
        String recombFile= dir+"13KRho.txt";
        newHomoSegDonorPrep(fileName,recombFile, .05, .3);
    }
}
