/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.FileInputStream;
import java.util.*;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.prefs.TasselPrefs;

/**
 *
 * @author kelly
 */
public class HMMForWholeBuild {
    public static String dir;
    public static byte diploidN= NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
    private static double[][] trans;//transition matrix, based on arbitrary f in constructor
    private static double[][] emiss;//emission matrix, based on the probability of observing a state given the read depth and HW expected frequency for the site
    private static MutableNucleotideAlignmentHDF5 mnah5;
    private static int[] whichWindow; //length of getSiteCount, hold window
    private static HashMap<Integer,Double> rho; //window->rho
    private static int[] windowOffset; //window->chromosome
    private static String[] whichChr;
    private static int[] sitesFromStartForward;//hold the last site since HMM resets to use for consensus HMM
    private static int[] sitesFromStartBackward;//hold the last site since HMM resets to use for consensus HMM
    private static boolean forward= true;
    private static byte[] mj;
    private static byte[] mn;
    private static double[] mnFreq;
    //*Takes in an hdf5 file, pulls out the heterozygous taxa, extract a matrix with depth info,*//
    //*calculate probabilities that actually homo/het for remaining sites based on HW expected *//
    //*frequencies and read depth. Does not take LD into account (assume independance for each site with respect to homo/het). *//
    //*Calculate joint probabilities based on 2pq for each site. Parse out the homozygous segments,*//
    //*returns the inbreds and inbred segments*//
    
    public static void run(String inFileRoot, String recombFile, String optionalTxtRoot, double hetCutoff, double f, boolean useTxtFile, boolean forwardAndReverse) {
        mnah5= MutableNucleotideAlignmentHDF5.getInstance(dir+inFileRoot+".hmp.h5");
        System.out.println("read in h5 file");
        String outFile= dir+inFileRoot+(useTxtFile==true?"InbredFor"+optionalTxtRoot:"InbredForHet"+hetCutoff)+".hmp.h5";
        ExportUtils.writeToMutableHDF5(mnah5, outFile, null, false);
        MutableNucleotideAlignmentHDF5 modMna= MutableNucleotideAlignmentHDF5.getInstance(outFile);
        System.out.println("generate output h5: "+outFile);
        readRecombinationFromTabTxt(recombFile, 10);
        if (mnah5.getLocus(0).getChromosomeName().equalsIgnoreCase(whichChr[0])==false||mnah5.getLocus(mnah5.getSiteCount()-1).getChromosomeName().equalsIgnoreCase(whichChr[whichChr.length-1])==false) 
            System.out.println("All chromosomes do not match between hapmap and recombination file. Sites with missing rho will be set to missing");
        MutableNucleotideAlignmentHDF5 mna;
        if (useTxtFile==true) {
            ParseAlignmentKelly.subsetHDF5FromTxt(dir, inFileRoot, optionalTxtRoot, false, false, false);
            System.out.println("subset alignment by text file");
            mna= MutableNucleotideAlignmentHDF5.getInstance(dir+inFileRoot+"StrictSubsetBy"+optionalTxtRoot+".hmp.h5");
        }
        else mna= MutableNucleotideAlignmentHDF5.getInstance(dir+inFileRoot+".hmp.h5");
        mj= new byte[mna.getSiteCount()];
        mn= new byte[mna.getSiteCount()];
        mnFreq= new double[mna.getSiteCount()];
        IdGroup subsetID= mna.getIdGroup();
        for (int site = 0; site < mna.getSiteCount(); site++) {
            mn[site]= mna.getMinorAllele(site);
            mj[site]= mna.getMajorAllele(site);
            mnFreq[site]= mna.getMinorAlleleFrequency(site);
        }
        System.out.println("generate mna");
        for (int taxon = 0; taxon < mnah5.getSequenceCount(); taxon++) {
            if (useTxtFile==false&&((double)mnah5.getHeterozygousCountForTaxon(taxon)/(double)mnah5.getTotalNotMissingForTaxon(taxon))<hetCutoff) continue;//if an inbred line, don't do anything
            if (useTxtFile==true&&(subsetID.whichIdNumber(mnah5.getFullTaxaName(taxon))<0)) continue;
            System.out.println("Working on taxon "+mnah5.getTaxaName(taxon));
            forward= true;
            sitesFromStartForward= new int[mnah5.getSiteCount()];
            byte[] currForward= newHomoSegDonorPrep(taxon,f);
            if (forwardAndReverse==true) {
                forward= false;
                sitesFromStartBackward= new int[mnah5.getSiteCount()];
                byte[] currReverse= newHomoSegDonorPrep(taxon,f);
                for (int site = 0; site < currForward.length; site++) { //if forward/reverse differ set to the longer path length
                    if (AlignmentUtils.isEqual(currForward[site],currReverse[site])) continue;
                    if (sitesFromStartBackward[site]>sitesFromStartForward[site]) currForward[site]= currReverse[site];
                }
            }
            modMna.setAllBases(taxon, currForward);
        }
        modMna.clean();
        
    }
    public static byte[] newHomoSegDonorPrep(int taxon, double f) {//this one only acts on LD and readDepth
        byte[] currAlign= mnah5.getBaseRange(taxon, 0, mnah5.getSiteCount());
        double probTrueHomo= -1;
        double probTrueHet= -1;
        boolean newChr= true;
        int lastSite= 0;
        int path= 0;
        for (int s = 0; s < mnah5.getSiteCount(); s++) {
            int site= forward==true?s:mnah5.getSiteCount()-1-s; //so the sites are correct for forward and reverse
            if (Arrays.binarySearch(whichChr, mnah5.getLocusName(site))<0) {//set to missing sites on chromosomes with no recombination information
                currAlign[site]= diploidN;
                continue;
            }
            if (mnFreq[site]==0) {//skip invariant sites and set to null if last site was called het and within 100 bp
                if (probTrueHet>probTrueHomo&&(Math.abs(mnah5.getSiteOfPhysicalPosition(site, mnah5.getLocus(site))
                        -mnah5.getSiteOfPhysicalPosition(lastSite, mnah5.getLocus(lastSite))))<100) currAlign[site]= diploidN;
                continue;
            } 
            if (mnah5.getBase(taxon, site)==diploidN) continue; //skip sites with no data
            if (mnah5.getLocus(site).getChromosomeName().equalsIgnoreCase(mnah5.getLocus(lastSite).getChromosomeName())==false) newChr= true;//reinitialize hmm for new chromosome
            int obs= mnah5.isHeterozygous(taxon, site)==true?1:0;
            double currRD= getReadDepthForMinMaj(mnah5.getDepthForAlleles(taxon, site),
                    getIndexForMinMaj(mnah5.getMajorAllele(site),mnah5.getMinorAllele(site)));
            double currPQ= 2*mj[site]*mn[site];
            //initialize emission matrix using HW/RD joint probability. The probability of observing a homo/het given both haplotypes are ibd {.999,.001} 
            //(sequencing error). The probability of observing a homo/het given 2 haplotypes present and these haplotypes vary at this particular site and rd>1
            //{(.5^read depth * 2pq)/2pq, (1-.5^read depth) * 2pq)/2pq}. Accounts for sampling based on read depth and assumes HW proportions for frequencies
            emiss= new double[][] {
                {.9999,.0001},//homo {obsHomo,obsHet}
                {currRD>1?(Math.pow(.5,currRD)*currPQ)/(currPQ):.99999,
                    currRD>1?((1-Math.pow(.5,currRD))*currPQ)/(currPQ):0.00001}//het {obsHomo,obsHet}
            };
            //the probability of recombination. could use the joint probability of recombination and inbreeding (f is made up) but too much bias to het (if use f=.2 from hapmap2)
            int dist= Math.abs(mnah5.getPositionInLocus(site)-mnah5.getPositionInLocus(lastSite));
            double currRho= (rho.get(whichWindow[site]));
            //prob of no recombination based on rho (pop estimate from jeff r-i hapmap2; crossovers/gen/bp). Assume one large population with whatever coalescent jeff assumed
            //double noRecomb= (dist>50)?(Math.pow((1-currRho),dist))+(IntMath.binomial(dist, 2)*Math.pow(currRho, 2)*(Math.pow((1-currRho),dist-2)))+(IntMath.binomial(dist, 4)*Math.pow(currRho, 4)*(Math.pow((1-currRho),dist-4))):Math.pow((1-currRho),dist);
            double noRecomb= Math.pow((1-currRho), dist);
            if (noRecomb<.0000000001) {//reset if sites are very far apart
                newChr= true;
                path= 0;
            } 
            if (newChr==true) path= 0;
            if (forward==true) sitesFromStartForward[site]= path;
            else sitesFromStartBackward[site]= path;
            
            if (obs==1) {//if obs heterozygous prob set to sequence error
                probTrueHomo= 0.001;
                probTrueHet= .999;
                lastSite= site;
                newChr= false;
            }
            else if (newChr==true) { //initialize HMM with HW expected allele freq (do if sites too far apart as well)
                probTrueHomo= (1-currPQ)*emiss[0][obs]; //the starting probabilities for each true state, given the observed state
                probTrueHet= currPQ*emiss[1][obs];
                lastSite= site;
                newChr= false;
            }
            else {
                trans= new double[][] {//probability of recombining and probability that recombination will bring together ibd haplotypes are assumed independant
                    {noRecomb,(1-noRecomb)*(1-f)},//homo {homo,het}
                    {(1-noRecomb)*f,noRecomb}//het {homo,het}
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
            if (probTrueHet>probTrueHomo) currAlign[site]= diploidN;
            path++;
        }
        return currAlign;
    }
    
    private static byte[] reverseAlignment(int taxon) { //reverses all the sites for a given taxon
        byte[] align= mnah5.getBaseRange(taxon, 0, mnah5.getSiteCount());
        byte[] reverse= new byte[mnah5.getSiteCount()];
        for (int site = 0; site < mnah5.getSiteCount(); site++) {
            reverse[site]= align[mnah5.getSiteCount()-site-1];
        }
        return reverse;
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
        TasselPrefs.putAlignmentRetainRareAlleles(false);
//        dir= "/Users/kls283/Desktop/";
//        dir= "/Users/kelly/Documents/GBS/Imputation/"; //laptopOld
        dir= "/home/kls283/Documents/Imputation/"; //desktop
        String fileNameRoot= "AllZeaGBS_v2.7wDepth";
        String recombFile= dir+"13KRho.txt";
        String txtRoot= "12S_RIMMA_Span_SEED";
        
//        String[] part= new String[]{"11","14","15","16","17"};
//        for (String i:part) {
//            run("AllZeaGBS_v2.7_SeqToGenos_part"+i+".hmp.h5", txt, recombFile, .005, .2, false, true);
//            run("AllZeaGBS_v2.7_SeqToGenos_part"+i+".hmp.h5", txt, recombFile, .005, .2, false, false);
//        }
        run(fileNameRoot,recombFile, txtRoot, .001, .2, true, false);
    }
}

