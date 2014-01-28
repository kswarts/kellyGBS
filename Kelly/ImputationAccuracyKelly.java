/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Scanner;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGenerator;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.prefs.TasselPrefs;
import org.apache.commons.lang.ArrayUtils;

/**
 *
 * @author kls283
 */
public class ImputationAccuracyKelly {
    public static String dir;
    public static boolean[][] knownHetMask;
    public static boolean[][] calledHomoMajMask;
    public static boolean[][] calledHomoMinMask;
    public static boolean[][] calledMissingMask;
    public static boolean[] highCovTaxa;
    public static boolean[] highCovHets;
    public static boolean[] highCovInbreds;
    public static boolean[] outbred;
    public static boolean[] inbred;
    public static byte diploidN= NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
    private static byte knownBase;
    private static byte impBase;
    private static byte knownMaj;
    private static byte knownMin;
    private static boolean knownHet;
    private static byte[] impArray;
    private static byte[] knownArray;
    private static int[] matchTaxon;
    private static double[][] perSiteTaxon;
    private static int knownIndex;
    private static int knownSiteCount;
    
    
    public static void maskFileSample(String inFile, int sampleIntensity) {
        Alignment a= ImportUtils.readGuessFormat(inFile,false);
        MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
        MutableNucleotideAlignment mnaKey= MutableNucleotideAlignment.getInstance(a);
        for (int taxon = 0; taxon < mnaKey.getSequenceCount(); taxon++) {
            for (int site = 0; site < mnaKey.getSiteCount(); site++) {mnaKey.setBase(taxon, site, diploidN);}
        }
        for (int taxon= 0;taxon<a.getSequenceCount();taxon++) {            
            for (int site= taxon;site<a.getSiteCount();site+= sampleIntensity+taxon) {
                mna.setBase(taxon, site, diploidN);
                mnaKey.setBase(taxon, site, a.getBase(taxon, site));
            }
        }
        mna.clean();
        ExportUtils.writeToMutableHDF5(mna, inFile.substring(0, inFile.indexOf(".hmp"))+"_masked.hmp.h5");
        ExportUtils.writeToMutableHDF5(mnaKey, inFile.substring(0, inFile.indexOf(".hmp"))+"_maskKey.hmp.h5");
    }
    
    public static void maskFile55k(String maskFile, boolean gzMask, String knownFile, boolean gzKnown) {
        Alignment a= ImportUtils.readFromHapmap(dir+maskFile+(gzMask==true?".hmp.txt.gz":".hmp.txt"), null);
        MutableNucleotideAlignment mask= MutableNucleotideAlignment.getInstance(a);
        Alignment known= ImportUtils.readFromHapmap(dir+knownFile+(gzKnown==true?".hmp.txt.gz":".hmp.txt"), null);
        int[] knownPos= known.getPhysicalPositions();
        for (int site = 0; site < mask.getSiteCount(); site++) {
            int pos= mask.getPositionInLocus(site);
            if (Arrays.binarySearch(knownPos, pos)<0) continue;
            for (int taxon = 0; taxon < mask.getSequenceCount(); taxon++) {
                mask.setBase(taxon, site, diploidN);
            }
        }
        mask.clean();
        ExportUtils.writeToHapmap(mask, true, dir+maskFile+"_masked55k.hmp.txt.gz", '\t', null);
    }
    //for depths 4 or more, requires hets to be called by more than one for less depth allele
    public static void maskFileByDepth(String depthFile, String inFile, int depthToMask, int maskDenom, boolean h5, boolean exportDepth) {
        System.out.println("Depth file: "+depthFile);
        System.out.println("File to mask: "+inFile);
        System.out.println("Site depth to mask: "+depthToMask);
        System.out.println("Divisor for physical positions to be masked: "+maskDenom);
        Alignment a= ImportUtils.readGuessFormat(inFile, true);
        MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
        MutableNucleotideAlignment siteMna= MutableNucleotideAlignment.getInstance(a);
        System.out.println("Read in file to mask");
        MutableNucleotideAlignmentHDF5 mnah5= MutableNucleotideAlignmentHDF5.getInstance(depthFile);
        System.out.println("Read in depth file");
        int cnt= 0;
        String[] h5Taxa= new String[mnah5.getSequenceCount()];
        for (int taxon = 0; taxon < h5Taxa.length; taxon++) { h5Taxa[taxon]= mnah5.getFullTaxaName(taxon);}
        ArrayList<Identifier> remove= new ArrayList<>();
        for (int taxon = 0; taxon < mna.getSequenceCount(); taxon++) {
            int taxaCnt= 0;
            int h5Taxon= Arrays.binarySearch(h5Taxa, mna.getFullTaxaName(taxon));
            if (h5Taxon<0) {
                System.out.println("Problem matching taxon "+mna.getFullTaxaName(taxon)+". Not masked and excluded.");
                remove.add(mna.getIdGroup().getIdentifier(taxon));
                continue;
            }
            for (int site = 0; site < mna.getSiteCount(); site++) {
                siteMna.setBase(taxon, site, diploidN);
                int h5Site= mnah5.getSiteOfPhysicalPosition(mna.getPositionInLocus(site), mna.getLocus(site));
                byte[] currDepth= mnah5.getDepthForAlleles(h5Taxon, h5Site);
                int[] currMinMaj= getIndexForMinMaj(mnah5.getMajorAllele(h5Site),mnah5.getMinorAllele(h5Site));
                if (getReadDepthForAlleles(currDepth,currMinMaj)==depthToMask) {
                    if ((depthToMask>3&&((getReadDepthForAlleles(currDepth,currMinMaj[0])==1)||(getReadDepthForAlleles(currDepth,currMinMaj[1])!=1)))) continue;
                    if (mnah5.getPositionInLocus(h5Site)%maskDenom==0) {
                        mna.setBase(taxon, site, diploidN);
                        siteMna.setBase(taxon, site, mnah5.getBase(h5Taxon, h5Site));
                        taxaCnt++;
                        cnt++;
                    }
                }   
            }
            System.out.println(taxaCnt+" sites masked for "+mna.getTaxaName(taxon));
        }
        System.out.println(cnt+" sites masked at a depth of "+depthToMask+" (site numbers that can be divided by "+maskDenom+")");
        System.out.println(remove.size()+" taxa not masked due to taxa name mismatch and excluded from output files");
        System.out.println(mna.getSequenceCount()-remove.size()+" total taxa output");
        
        IdGroup removeIDs= IdGenerator.createIdGroup(remove.size());
        for (int i = 0; i < remove.size(); i++) {removeIDs.setIdentifier(i, remove.get(i));}
        mna.clean();
        siteMna.clean();
        Alignment filterMna= FilterAlignment.getInstanceRemoveIDs(mna, removeIDs);
        Alignment filterSiteMna= FilterAlignment.getInstanceRemoveIDs(siteMna, removeIDs);
        
        if (h5==true) {
            if (exportDepth==true) ExportUtils.writeToMutableHDF5((Alignment)filterMna, inFile.length()-6+"_maskedDepth"+depthToMask+"_Denom"+maskDenom, null, true);
            else {
                ExportUtils.writeToMutableHDF5(filterMna,inFile.substring(0, inFile.length()-6)+"_masked_Depth"+depthToMask+"_Denom"+maskDenom+".hmp.h5", null, false);
                ExportUtils.writeToMutableHDF5(filterSiteMna,inFile.substring(0, inFile.length()-6)+"_maskKey_Depth"+depthToMask+"_Denom"+maskDenom+".hmp.h5", null, false);
            }
        }
        else {
            ExportUtils.writeToHapmap(filterMna, true, inFile.substring(0, inFile.length()-6)+"_masked_Depth"+depthToMask+"_Denom"+maskDenom+".hmp.txt.gz", '\t', null);
            ExportUtils.writeToHapmap(filterSiteMna, true, inFile.substring(0, inFile.length()-6)+"_maskKey_Depth"+depthToMask+"_Denom"+maskDenom+".hmp.txt.gz", '\t', null);
        }
    }
    
    //for depths 4 or more, requires hets to be called by more than one for less depth allele
    public static void maskBaseFileByDepth(String depthFile, int depthToMask, int maskDenom, boolean exportDepth) {
        System.out.println("Depth file: "+depthFile);
        System.out.println("Site depth to mask: "+depthToMask);
        System.out.println("Divisor for physical positions to be masked: "+maskDenom);
        MutableNucleotideAlignmentHDF5 mnah5= MutableNucleotideAlignmentHDF5.getInstance(depthFile);
        System.out.println("Generate file to mask");
        String outMasked= depthFile.substring(0, depthFile.length()-7)+"_masked_Depth"+depthToMask+"_Denom"+maskDenom+".hmp.h5";
        String outKey= depthFile.substring(0, depthFile.length()-7)+"_maskKey_Depth"+depthToMask+"_Denom"+maskDenom+".hmp.h5";
        ExportUtils.writeToMutableHDF5(mnah5,outMasked, null, exportDepth);
        ExportUtils.writeToMutableHDF5(mnah5,outKey, null, false);
        MutableNucleotideAlignmentHDF5 mna= MutableNucleotideAlignmentHDF5.getInstance(outMasked);
        System.out.println("read back in maskedFile: "+outMasked);
        MutableNucleotideAlignmentHDF5 key= MutableNucleotideAlignmentHDF5.getInstance(outKey);
        System.out.println("read back in keyFile: "+outKey);
        int cnt= 0;
        System.out.println("Starting mask...");
        for (int taxon = 0; taxon < mna.getSequenceCount(); taxon++) {
            int taxaCnt= 0;
            byte[] taxonMask= new byte[mna.getSiteCount()];
            byte[] taxonKey= new byte[mna.getSiteCount()];
            for (int site = 0; site < mna.getSiteCount(); site++) {
                taxonKey[site]= diploidN;
                taxonMask[site]= mnah5.getBase(taxon, site);
                byte[] currDepth= mnah5.getDepthForAlleles(taxon, site);
                int[] currMinMaj= getIndexForMinMaj(mnah5.getMajorAllele(site),mnah5.getMinorAllele(site));
                if (getReadDepthForAlleles(currDepth,currMinMaj)==depthToMask) {
                    if ((AlignmentUtils.isHeterozygous(taxonMask[site])==true)&&(depthToMask>3)&&((getReadDepthForAlleles(currDepth,currMinMaj[0])==1)||(getReadDepthForAlleles(currDepth,currMinMaj[1])==1))) continue;
                    if (mnah5.getPositionInLocus(site)%maskDenom==0) {
                        taxonMask[site]= diploidN;
                        taxonKey[site]= mnah5.getBase(taxon, site);
                        taxaCnt++;
                        cnt++;
                    }
                }   
            }
            mna.setAllBases(taxon, taxonMask);
            key.setAllBases(taxon, taxonKey);
            
            System.out.println(taxaCnt+" sites masked for "+mna.getTaxaName(taxon));
        }
        System.out.println(cnt+" sites masked at a depth of "+depthToMask+" (site numbers that can be divided by "+maskDenom+")");
        mna.clean();
        key.clean();
    }
    
    //returns an array the length of key.getSiteCount() that contains the MAF class index of the donor file or -1 if not included
    public static int[] readInMAFFile (String donorMAFFile, Alignment key, double[] MAFClass) {
        int[] MAF= new int[key.getSiteCount()];
        for (int i = 0; i < MAF.length; i++) {MAF[i]= -1;} //initialize to -1
        try {
            FileInputStream fis= new FileInputStream(donorMAFFile);
            Scanner scanner= new Scanner(fis);
            do {
                String next= scanner.nextLine();
                if (next.isEmpty()) break;
                String[] vals= next.split("\t");
                int site= key.getSiteOfPhysicalPosition(Integer.parseInt(vals[1]),key.getLocus(vals[0]));
                if (site<0) continue;
                int currClass= Arrays.binarySearch(MAFClass, Double.parseDouble(vals[2]));
                MAF[site]= currClass<0?Math.abs(currClass)-1:currClass;
            }
            while (scanner.hasNextLine());
            scanner.close();
            fis.close();
        }
        catch (Exception e) {
            System.out.println("Problem reading in taxa names");
        }
        return MAF;
    }
    
    public static void imputeUsingPQ(String unimputedFileName) {
        Alignment orig = ImportUtils.readGuessFormat(unimputedFileName, false);
        String newFile= unimputedFileName.substring(0, unimputedFileName.indexOf(".hmp"))+"ImputedPQ.hmp.h5";
        ExportUtils.writeToMutableHDF5(orig, newFile);
        MutableNucleotideAlignmentHDF5 impute= MutableNucleotideAlignmentHDF5.getInstance(newFile);
        byte[] maj= new byte[impute.getSiteCount()];
        byte[] min= new byte[impute.getSiteCount()];
        double[] freqQ= new double[impute.getSiteCount()];
        for (int site = 0; site < impute.getSiteCount(); site++) {
            maj[site]= impute.getMajorAllele(site);
            min[site]= impute.getMinorAllele(site);
            freqQ[site]= impute.getMinorAlleleFrequency(site);
        }
        for (int taxon = 0; taxon < impute.getSequenceCount(); taxon++) {
            byte[] currAlign= impute.getBaseRow(taxon);
            for (int site = 0; site < impute.getSiteCount(); site++) {
                if (currAlign[site]!=Alignment.UNKNOWN_DIPLOID_ALLELE) continue;
                double rand= Math.random();
                double q2= (freqQ[site])*(freqQ[site]);
                double twopq= 2*freqQ[site]*(1-freqQ[site]);
                if (rand<q2) currAlign[site]= AlignmentUtils.getDiploidValue(min[site], min[site]);
                else if (rand<(q2+twopq)) currAlign[site]= AlignmentUtils.getDiploidValue(maj[site], min[site]);
                else currAlign[site]= AlignmentUtils.getDiploidValue(maj[site], maj[site]);
            }
            impute.setAllBases(taxon, currAlign);
        }
        impute.clean();
    }
    
    //MAFClass should always have a non-polymorphic 0 class
     private static int[] localMAF (String unimputedFileName, Alignment imputed, double[] MAFClass) {
         Alignment unimputed= ImportUtils.readGuessFormat(unimputedFileName, true);
         int[] MAF= new int[imputed.getSiteCount()];
         for (int i = 0; i < MAF.length; i++) {MAF[i]= -1;} //initialize to -1
         for (int site = 0; site < imputed.getSiteCount(); site++) {
             int unimpSite= unimputed.getSiteOfPhysicalPosition(imputed.getPositionInLocus(site), imputed.getLocus(site));
             if (unimpSite<0) continue;
             double maf= unimputed.getMinorAlleleFrequency(unimpSite);
             int currClass= Arrays.binarySearch(MAFClass, maf);
             MAF[site]= currClass<0?Math.abs(currClass)-1:currClass;
         }
        return MAF;
    }
    
    public static void accuracyDepth(String keyFile, String imputedFileName, String unimputedFileName, String donorMAFFile, double[] MAFClass) {
        System.out.println("Key file: "+keyFile+"\nUnimputed file: "+unimputedFileName+"\nImputed file: "+imputedFileName);
        Alignment index = ImportUtils.readGuessFormat(keyFile, false);
        System.out.println("\nKey file read in");
        Alignment imputed= (imputedFileName.contains(".vcf"))?ImportUtils.readFromVCF(imputedFileName, null, 2):ImportUtils.readGuessFormat(imputedFileName, true);
        System.out.println("Imputed file read in");
        double[][] all= new double[3][5]; //arrays held ("columns"): 0-maskedMinor, 1-maskedHet, 2-maskedMajor; each array ("rows"):0-to minor, 1-to het, 2-to major, 3-unimp, 4-total for known type
        DecimalFormat df = new DecimalFormat("0.########");
        if (Arrays.equals(index.getPhysicalPositions(), imputed.getPhysicalPositions())==false) System.out.println("Not all physical positions match between imputed alignment and key. Only use those that match");
        boolean MAFon= (MAFClass==null||donorMAFFile==null)?false:true;
        int[] MAF= null; double[][][] mafAll= null;//the first array holds the separate MAF classes, and the second and third hold the results (like all)
        if (MAFon==true) {//If provided with a file to read in (from donor files)
            mafAll= new double[MAFClass.length][3][5];
            MAF= readInMAFFile(donorMAFFile, imputed, MAFClass);//MAF sites match the imputed taxon
            System.out.println("MAF taken from input txt file");
        }
        else if (unimputedFileName!= null&&MAFClass!=null) {//if it should read MAF from the unimputed file (ie out of beagle)
            MAF= localMAF(unimputedFileName,imputed,MAFClass);
            System.out.println("Unimputed file read in and MAF generated from allele frequencies of this file");
            MAFon= true; 
            mafAll= new double[MAFClass.length][3][5];
        }
        ArrayList<Double> x= new ArrayList<>();
        ArrayList<Double> y= new ArrayList<>();
        int keySite= 0;
        int keyTaxon= 0;
        IdGroup keyTaxaNames= index.getIdGroup();
        for (int taxon = 0; taxon < imputed.getSequenceCount(); taxon++) {
            String currTaxon= imputed.getFullTaxaName(taxon);
            keyTaxon= keyTaxaNames.whichIdNumber(currTaxon);
            if (keyTaxon<0) continue;
            boolean currMAF;
            int maf= -1; //holds the MAF class for the current site
            for (int site = 0; site < imputed.getSiteCount(); site++) {
                if (MAFon) maf= MAF[site];
                currMAF= (MAFon&&maf>-1)?true:false;
                keySite= index.getSiteOfPhysicalPosition(imputed.getPositionInLocus(site), index.getLocus(imputed.getLocusName(site)));
                if (keySite<0) continue;
                byte known = index.getBase(keyTaxon, keySite);
                String knownBase= index.getBaseAsString(keyTaxon, keySite);
                if (known == diploidN) continue;
                byte imp = imputed.getBase(taxon, site);
                String impBase= imputed.getBaseAsString(taxon,site);
                if (AlignmentUtils.isHeterozygous(known) == true) {
                    all[1][4]++;
                    if (currMAF) mafAll[maf][1][4]++;
                    if (imp == diploidN) {all[1][3]++; if (currMAF) mafAll[maf][1][3]++;}
                    else if (AlignmentUtils.isEqual(imp, known) == true) {all[1][1]++;x.add(1.0);y.add(1.0);if (currMAF) mafAll[maf][1][1]++;}
                    else if (AlignmentUtils.isHeterozygous(imp) == false && AlignmentUtils.isPartiallyEqual(imp, imputed.getMinorAllele(site)) == true) {//to minor
                        x.add(1.0);y.add(0.0);
                        all[1][0]++;
                        if (currMAF) mafAll[maf][1][0]++;
                    }
                    else if (AlignmentUtils.isHeterozygous(imp) == false && AlignmentUtils.isPartiallyEqual(imp, imputed.getMajorAllele(site)) == true){
                        all[1][2]++;
                        if (currMAF) mafAll[maf][1][2]++;
                        x.add(1.0);y.add(2.0);
                    }
                    else {System.out.println("Exclude: More than two allele states at site index "+site); all[1][4]--;if (currMAF) mafAll[maf][1][4]--;}
                } 
                else if (index.getBaseArray(keyTaxon, keySite)[0] == imputed.getMinorAllele(site)) {
                    all[0][4]++;
                    if (currMAF) mafAll[maf][0][4]++;
                    if (imp == diploidN) {all[0][3]++; if (currMAF) mafAll[maf][0][3]++;}
                    else if (AlignmentUtils.isEqual(imp, known) == true) {all[0][0]++;x.add(0.0);y.add(0.0);if (currMAF) mafAll[maf][0][0]++;}
                    else if (AlignmentUtils.isHeterozygous(imp) == true && AlignmentUtils.isPartiallyEqual(imp, known) == true) {
                        all[0][1]++;
                        x.add(0.0);y.add(1.0);
                        if (currMAF) mafAll[maf][0][1]++;
                    }
                    else {
                        all[0][2]++;
                        x.add(0.0);y.add(2.0);
                        if (currMAF) mafAll[maf][0][2]++;
                    }
                }
                else if (index.getBaseArray(keyTaxon, keySite)[0] == imputed.getMajorAllele(site)) {
                    all[2][4]++;
                    if (currMAF) mafAll[maf][2][4]++;
                    if (imp == diploidN) {all[2][3]++; if (currMAF) mafAll[maf][2][3]++;}
                    else if (AlignmentUtils.isEqual(imp, known) == true) {all[2][2]++;x.add(2.0);y.add(2.0);if (currMAF) mafAll[maf][2][2]++;}
                    else if (AlignmentUtils.isHeterozygous(imp) == true && AlignmentUtils.isPartiallyEqual(imp, known) == true) {all[2][1]++;x.add(2.0);y.add(1.0);if (currMAF) mafAll[maf][2][1]++;}
                    else {all[2][0]++;x.add(2.0);y.add(0.0);if (currMAF) mafAll[maf][2][0]++;}
                }
                else continue;
            }
        }
        double r2= pearsonR2(x,y);
        //for (int i = 0; i < x.size(); i++) {System.out.println(x.get(i)+"\t"+y.get(i));} //to output correlation xy
        double[][] forLambda= new double[3][3];
        forLambda[0]= ArrayUtils.subarray(all[0], 0, 3);forLambda[1]= ArrayUtils.subarray(all[1], 0, 3);forLambda[2]=ArrayUtils.subarray(all[2], 0, 3);
        double l= lambda(forLambda);
        try {
            File outputFile = new File(imputedFileName.substring(0, imputedFileName.length()-7) + "DepthAccuracy.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("##Taxon\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumMinor\tCorrectMinor\tMinorToHet\tMinorToMajor\tUnimpMinor"
                    + "\tNumHets\tHetToMinor\tCorrectHet\tHetToMajor\tUnimpHet\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimpMajor\tr2\tlambda\n");
            outStream.writeBytes("##TotalByImputed\t"+(all[0][4]+all[1][4]+all[2][4])+"\t"+(all[0][4]+all[1][4]+all[2][4]-all[0][3]-all[1][3]-all[2][3])+"\t"+
                    ((all[0][3]+all[1][3]+all[2][3])/(all[0][4]+all[1][4]+all[2][4]))+"\t"+all[0][4]+"\t"+all[0][0]+"\t"+all[0][1]+"\t"+all[0][2]+"\t"+all[0][3]+
                    "\t"+all[1][4]+"\t"+all[1][0]+"\t"+all[1][1]+"\t"+all[1][2]+"\t"+all[1][3]+"\t"+all[2][4]+"\t"+all[2][0]+"\t"+all[2][1]+"\t"+all[2][2]+
                    "\t"+all[2][3]+"\t"+r2+"\t"+l+"\n");
            outStream.writeBytes("#Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nx\ty\tN\tprop\n"
                    +0+"\t"+0+"\t"+all[0][0]+"\t"+df.format((all[0][0])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +0+"\t"+.5+"\t"+all[0][1]+"\t"+df.format((all[0][1])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +0+"\t"+1+"\t"+all[0][2]+"\t"+df.format((all[0][2])/(all[0][0]+all[0][1]+all[0][2]))+"\n"
                    +.5+"\t"+0+"\t"+all[1][0]+"\t"+df.format((all[1][0])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +.5+"\t"+.5+"\t"+all[1][1]+"\t"+df.format((all[1][1])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +.5+"\t"+1+"\t"+all[1][2]+"\t"+df.format((all[1][2])/(all[1][0]+all[1][1]+all[1][2]))+"\n"
                    +1+"\t"+0+"\t"+all[2][0]+"\t"+df.format((all[2][0])/(all[2][0]+all[2][1]+all[2][2]))+"\n"
                    +1+"\t"+.5+"\t"+all[2][1]+"\t"+df.format((all[2][1])/(all[2][0]+all[2][1]+all[2][2]))+"\n"
                    +1+"\t"+1+"\t"+all[2][2]+"\t"+df.format((all[2][2])/(all[2][0]+all[2][1]+all[2][2]))+"\n");
            outStream.writeBytes("#Proportion unimputed:\n#minor <- "+all[0][3]/all[0][4]+"\n#het<- "+all[1][3]/all[1][4]+"\n#major<- "+all[2][3]/all[2][4]);
            System.out.println("##Taxon\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumMinor\tCorrectMinor\tMinorToHet\tMinorToMajor\tUnimpMinor"
                    + "\tNumHets\tHetToMinor\tCorrectHet\tHetToMajor\tUnimpHet\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimpMajor\tr2\tlambda");
            System.out.println("TotalByImputed\t"+(all[0][4]+all[1][4]+all[2][4])+"\t"+(all[0][4]+all[1][4]+all[2][4]-all[0][3]-all[1][3]-all[2][3])+"\t"+
                    ((all[0][3]+all[1][3]+all[2][3])/(all[0][4]+all[1][4]+all[2][4]))+"\t"+all[0][4]+"\t"+all[0][0]+"\t"+all[0][1]+"\t"+all[0][2]+"\t"+all[0][3]+
                    "\t"+all[1][4]+"\t"+all[1][0]+"\t"+all[1][1]+"\t"+all[1][2]+"\t"+all[1][3]+"\t"+all[2][4]+"\t"+all[2][0]+"\t"+all[2][1]+"\t"+all[2][2]+
                    "\t"+all[2][3]+"\t"+r2+"\t"+l);
            System.out.println("Proportion unimputed:\nminor: "+all[0][3]/all[0][4]+"\nhet: "+all[1][3]/all[1][4]+"\nmajor: "+all[2][3]/all[2][4]);
            System.out.println("#Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nx\ty\tN\tprop\n"
                    +0+"\t"+0+"\t"+all[0][0]+"\t"+(all[0][0])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +0+"\t"+.5+"\t"+all[0][1]+"\t"+(all[0][1])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +0+"\t"+1+"\t"+all[0][2]+"\t"+(all[0][2])/(all[0][0]+all[0][1]+all[0][2])+"\n"
                    +.5+"\t"+0+"\t"+all[1][0]+"\t"+(all[1][0])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +.5+"\t"+.5+"\t"+all[1][1]+"\t"+(all[1][1])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +.5+"\t"+1+"\t"+all[1][2]+"\t"+(all[1][2])/(all[1][0]+all[1][1]+all[1][2])+"\n"
                    +1+"\t"+0+"\t"+all[2][0]+"\t"+(all[2][0])/(all[2][0]+all[2][1]+all[2][2])+"\n"
                    +1+"\t"+.5+"\t"+all[2][1]+"\t"+(all[2][1])/(all[2][0]+all[2][1]+all[2][2])+"\n"
                    +1+"\t"+1+"\t"+all[2][2]+"\t"+(all[2][2])/(all[2][0]+all[2][1]+all[2][2])+"\n");
            outStream.close();
        } catch (Exception e) {
            System.out.println(e);
        }
        if (MAFon) try {
            File outputFile = new File(imputedFileName.substring(0, imputedFileName.length()-7) + "DepthAccuracyMAF.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("##\tMAFClass\tTotalSitesMasked\tTotalSitesCompared\tTotalPropUnimputed\tNumHets\tHetToMinor\tHetToMajor\tCorrectHet\tUnimpHet\tNumMinor\tMinorToMajor\tMinorToHet\tCorrectMinor\t"
                    + "UnimpMinor\tNumMajor\tMajorToMinor\tMajorToHet\tCorrectMajor\tUnimputedMajor\tr2\tlambda\n");
            for (int i= 0; i<MAFClass.length;i++) {
                outStream.writeBytes("##TotalByImputed\t"+MAFClass[i]+"\t"+(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4])+"\t"+(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4]-mafAll[i][0][3]-mafAll[i][1][3]-mafAll[i][2][3])+"\t"+
                    ((mafAll[i][0][3]+mafAll[i][1][3]+mafAll[i][2][3])/(mafAll[i][0][4]+mafAll[i][1][4]+mafAll[i][2][4]))+"\t"+mafAll[i][0][4]+"\t"+mafAll[i][0][0]+"\t"+mafAll[i][0][1]+"\t"+mafAll[i][0][2]+"\t"+mafAll[i][0][3]+
                    "\t"+mafAll[i][1][4]+"\t"+mafAll[i][1][0]+"\t"+mafAll[i][1][1]+"\t"+mafAll[i][1][2]+"\t"+mafAll[i][1][3]+"\t"+mafAll[i][2][4]+"\t"+mafAll[i][2][0]+"\t"+mafAll[i][2][1]+"\t"+mafAll[i][2][2]+
                    "\t"+mafAll[i][2][3]+"\t"+r2+"\t"+l+"\n");
            }
            outStream.writeBytes("#MAFClass,Minor=0,Het=1,Major=2;x is masked(known), y is predicted\nMAF\tx\ty\tN\tprop\n");
            for (int i= 0; i<MAFClass.length;i++) { outStream.writeBytes(
                    MAFClass[i]+"\t"+0+"\t"+0+"\t"+mafAll[i][0][0]+"\t"+df.format((mafAll[i][0][0])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+0+"\t"+.5+"\t"+mafAll[i][0][1]+"\t"+df.format((mafAll[i][0][1])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+0+"\t"+1+"\t"+mafAll[i][0][2]+"\t"+df.format((mafAll[i][0][2])/(mafAll[i][0][0]+mafAll[i][0][1]+mafAll[i][0][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+0+"\t"+mafAll[i][1][0]+"\t"+df.format((mafAll[i][1][0])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+.5+"\t"+mafAll[i][1][1]+"\t"+df.format((mafAll[i][1][1])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+.5+"\t"+1+"\t"+mafAll[i][1][2]+"\t"+df.format((mafAll[i][1][2])/(mafAll[i][1][0]+mafAll[i][1][1]+mafAll[i][1][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+0+"\t"+mafAll[i][2][0]+"\t"+df.format((mafAll[i][2][0])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+.5+"\t"+mafAll[i][2][1]+"\t"+df.format((mafAll[i][2][1])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n"
                    +MAFClass[i]+"\t"+1+"\t"+1+"\t"+mafAll[i][2][2]+"\t"+df.format((mafAll[i][2][2])/(mafAll[i][2][0]+mafAll[i][2][1]+mafAll[i][2][2]))+"\n");
            }
            outStream.writeBytes("#Proportion unimputed:\n#MAF\tminor\thet\tmajor\n");
            for (int i= 0; i<MAFClass.length;i++) { 
                outStream.writeBytes("#"+MAFClass[i]+"\t"+mafAll[i][0][3]/mafAll[i][0][4]+"\t"+mafAll[i][1][3]/mafAll[i][1][4]+"\t"+mafAll[i][2][3]/mafAll[i][2][4]+"\n");
            }
            outStream.flush();
            outStream.close();
        } catch (Exception e) {
            System.out.println(e);
        }
    }
    
    public static void jointAccuracyDepth(String keyFile, String imputedFileName, String imputedBeagleFile) {
        Alignment index = ImportUtils.readGuessFormat(keyFile, false);
        Alignment imputed = ImportUtils.readGuessFormat(imputedFileName, false);
        Alignment beagle= (imputedBeagleFile.contains(".vcf"))?ImportUtils.readFromVCF(imputedBeagleFile, null, 6):ImportUtils.readGuessFormat(imputedBeagleFile, false);
        String outName= imputedBeagleFile.substring(0,imputedBeagleFile.indexOf("BeagleBase"))+"Consensus.hmp.h5";
        ExportUtils.writeToMutableHDF5(imputed, outName);
        MutableNucleotideAlignmentHDF5 con= MutableNucleotideAlignmentHDF5.getInstance(outName);
        double[][] all= new double[3][13]; //arrays held ("columns"): 0-maskedMinor, 1-maskedHet, 2-maskedMajor; each array ("rows" Beagle/Ours):
        //0-minor/minor, 1-minor/het, 2-minor/major, 3-minor/unimp, 4-het/minor, 5-het/het, 6-het/major, 7-het/unimputed, 8-major/minor, 9-major/het, 
        //10-major/major, 11-major/unimputed, 12-total for known type
        String[] classes= new String[]{"Minor/Minor","Minor/Het","Minor/Major","Minor/Unimputed","Het/Minor","Het/Het","Het/Major","Het/Unimputed",
        "Major/Minor","Major/Het","Major/Major","Major/Unimputed","TotalKnown"};
        DecimalFormat df = new DecimalFormat("0.########");
        System.out.println("Imputed file: "+imputedFileName+"\nKey file: "+keyFile);
        if (Arrays.equals(index.getPhysicalPositions(), imputed.getPhysicalPositions())==false) System.out.println("sites do not match between imputed alignment and key");
        int keySite= 0; int keyTaxon= 0; int beagleSite= 0; int beagleTaxon= 0;
        IdGroup keyTaxaNames= index.getIdGroup();
        byte[] empty= new byte[imputed.getSiteCount()];
        for (int i = 0; i < empty.length; i++) {empty[i]= Alignment.UNKNOWN_DIPLOID_ALLELE;}
        for (int taxon = 0; taxon < imputed.getSequenceCount(); taxon++) {
            byte[] newSeq= empty;
            String currTaxon= imputed.getFullTaxaName(taxon);
            keyTaxon= keyTaxaNames.whichIdNumber(currTaxon);
            beagleTaxon= beagle.getIdGroup().whichIdNumber(currTaxon);
            if (beagleTaxon<0) continue;
            if (keyTaxon<0) continue;
            for (int site = 0; site < imputed.getSiteCount(); site++) {
                keySite= index.getSiteOfPhysicalPosition(imputed.getPositionInLocus(site), imputed.getLocus(site));
                byte known = index.getBase(keyTaxon, keySite);
                String knownBase= index.getBaseAsString(keyTaxon, keySite);
                if (known == diploidN) continue;
                byte maj= imputed.getMajorAllele(site);
                byte imp = imputed.getBase(taxon, site);
                String impBase= imputed.getBaseAsString(taxon,site);
                int impLocus= imputed.getPositionInLocus(site);
                beagleSite= beagle.getSiteOfPhysicalPosition(imputed.getPositionInLocus(site), beagle.getLocus(imputed.getLocusName(site)));
                if (beagleSite<0) continue;
                byte beagleBase= beagle.getBase(beagleTaxon, beagleSite);
                String beagleBaseString= beagle.getBaseAsString(beagleTaxon, beagleSite);
                if (AlignmentUtils.isHeterozygous(known) == true) {//known is het
                    all[1][12]++;
                    if (AlignmentUtils.isHeterozygous(beagleBase)==true) {//beagle is het
                        if (AlignmentUtils.isHeterozygous(imp)==true) {all[1][5]++; newSeq[site]= beagleBase;}
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[1][7]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) all[1][6]++;
                        else all[1][4]++;
                    }
                    else if (AlignmentUtils.isPartiallyEqual(beagleBase, maj)==false) {//beagle is minor
                        if (AlignmentUtils.isHeterozygous(imp)==true) all[1][1]++;
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[1][3]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) all[1][2]++;
                        else {all[1][0]++; newSeq[site]= beagleBase;}
                    }
                    else if (AlignmentUtils.isPartiallyEqual(beagleBase, maj)==true) {//beagle is major
                        if (AlignmentUtils.isHeterozygous(imp)==true) all[1][9]++;
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[1][11]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) {all[1][10]++; newSeq[site]= beagleBase;}
                        else all[1][8]++;
                    }
                    else {System.out.println("Exclude: More than two allele states at position "+imputed.getPositionInLocus(site)); all[1][12]--;}
                }
                else if (AlignmentUtils.isPartiallyEqual(known, maj) == true) {//known is major
                    all[2][12]++;
                    if (AlignmentUtils.isHeterozygous(beagleBase)==true) {//beagle is het
                        if (AlignmentUtils.isHeterozygous(imp)==true) {all[2][5]++; newSeq[site]= beagleBase;}
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[2][7]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) all[2][6]++;
                        else all[2][4]++;
                    }
                    else if (AlignmentUtils.isPartiallyEqual(beagleBase, maj)==false) {//beagle is minor
                        if (AlignmentUtils.isHeterozygous(imp)==true) all[2][1]++;
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[2][3]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) all[2][2]++;
                        else {all[2][0]++; newSeq[site]= beagleBase;}
                    }
                    else if (AlignmentUtils.isPartiallyEqual(beagleBase, maj)==true) {//beagle is major
                        if (AlignmentUtils.isHeterozygous(imp)==true) all[2][9]++;
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[2][11]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) {all[2][10]++; newSeq[site]= beagleBase;}
                        else all[2][8]++;
                    }
                    else {System.out.println("Exclude: More than two allele states at position "+imputed.getPositionInLocus(site)); all[2][12]--;}
                }
                else if (AlignmentUtils.isPartiallyEqual(known, maj) == false) {//known is minor
                    all[0][12]++;
                    if (AlignmentUtils.isHeterozygous(beagleBase)==true) {//beagle is het
                        if (AlignmentUtils.isHeterozygous(imp)==true) {all[0][5]++; newSeq[site]= beagleBase;}
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[0][7]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) all[0][6]++;
                        else all[0][4]++;
                    }
                    else if (AlignmentUtils.isPartiallyEqual(beagleBase, maj)==false) {//beagle is minor
                        if (AlignmentUtils.isHeterozygous(imp)==true) all[0][1]++;
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[0][3]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) all[0][2]++;
                        else {all[0][0]++; newSeq[site]= beagleBase;}
                    }
                    else if (AlignmentUtils.isPartiallyEqual(beagleBase, maj)==true) {//beagle is major
                        if (AlignmentUtils.isHeterozygous(imp)==true) all[0][9]++;
                        else if (AlignmentUtils.isEqual(imp, Alignment.UNKNOWN_DIPLOID_ALLELE)) all[0][11]++;
                        else if (AlignmentUtils.isPartiallyEqual(imp, maj)) {all[0][10]++; newSeq[site]= beagleBase;}
                        else all[0][8]++;
                    }
                    else {System.out.println("Exclude: More than two allele states at position "+imputed.getPositionInLocus(site)); all[0][12]--;}
                }
                else continue;
            }
            con.setBaseRange(taxon, 0, newSeq);
        }
        System.out.println("\tKnownMinor\tKnownHet\tKnownMajor");
        for (int i = 0; i < all[0].length; i++) {System.out.println(classes[i]+"\t"+all[0][i]+"\t"+all[1][i]+"\t"+all[2][i]);}
        try {
            File outputFile = new File(imputedFileName.substring(0, imputedFileName.length()-7) + "JointAccuracy.txt");
            DataOutputStream outStream = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFile)));
            outStream.writeBytes("\tKnownMinor\tKnownHet\tKnownMajor");
            for (int i = 0; i < all[0].length; i++) {outStream.writeBytes("\n"+classes[i]+"\t"+all[0][i]+"\t"+all[1][i]+"\t"+all[2][i]);}
            outStream.close();
        } catch (Exception e) {
            System.out.println(e);
        }
        con.clean();
    }
    
    //this is the sample multiple r2. input should be coded if categorical according to desired linear contrast
    private static double pearsonR2(ArrayList<Double> x, ArrayList<Double> y) {
        double[][] xy= new double[2][x.size()];
        xy[0]= ArrayUtils.toPrimitive(x.toArray(new Double[x.size()]), -1.0); xy[1]= ArrayUtils.toPrimitive(y.toArray(new Double[y.size()]), -1.0);
        double meanX= 0; double meanY= 0; double varX= 0; double varY= 0; double covXY= 0; double r2= 0.0;
        for (int i = 0; i < xy[0].length; i++) {meanX+=xy[0][i]; meanY+= xy[1][i];}
        meanX= meanX/(xy[0].length-1); meanY= meanY/(xy[1].length-1);
        double currX, currY;
        for (int i = 0; i < xy[0].length; i++) {
            currX= xy[0][i]-meanX; currY= xy[1][i]-meanY;
            varX+= currX*currX; varY+= currY*currY;
            covXY+= currX*currY;
        }
        r2= (covXY/(Math.sqrt(varX)*Math.sqrt(varY)))*(covXY/(Math.sqrt(varX)*Math.sqrt(varY)));
        return r2;
    }
    
    //implements lambda, a non-parametric test of association for categorical nominal variables
    //input: first array is independant, second is dependant (double[2][num classes])
    private static double lambda(double[][] id) {
        double[] totalDep= new double[id[0].length+1];
        double e1= 0;
        for (int j = 0; j < id.length; j++) {
            for (int i = 0; i < id[0].length; i++) { totalDep[j]+= id[i][j];}
            if (totalDep[j]>e1) e1= totalDep[j];
        }
        double e2= 0;
        for (int i = 0; i < id.length; i++) {
            double currMode= 0;
            for (int j = 0; j < id[0].length; j++) {
                if (id[i][j]>currMode) currMode= id[i][j];
            }
            e2+= currMode;
        }
        double l= (e1-e2)/e1;
        return l;
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
    
    public static int getReadDepthForAlleles(byte[] depthForAlleles, int[] alleleIndices) {
        int depth= 0;
        for (int all: alleleIndices) {
            depth+= (int)depthForAlleles[all];
        }
        return depth;
    }
    
    public static int getReadDepthForAlleles(byte[] depthForAlleles, int allele) {
        int depth= (int)depthForAlleles[allele];
        return depth;
    }
    
    public static double getTrue(boolean[] mask) {
        int count= 0;
        for (int i= 0;i<mask.length;i++) {
            if (mask[i]==true) count++;
        }
        return count;
    }  
    
    public static double getTrue(boolean[][] mask) {
        int count= 0;
        for (int i= 0;i<mask.length;i++) {
            for (int j= 0;j<mask[0].length;j++) {
                if (mask[i][j]==true) count++;
            }   
        }
        return count;
    }
    
    public static double getTrue(boolean[][] mask1, boolean[] mask2) {
        int count= 0;
        for (int i= 0;i<mask1.length;i++) {
            for (int j= 0;j<mask1[0].length;j++) {
                if (mask1[i][j]==true&&mask2[i]==true) count++;
            }   
        }
        return count;
    }
    
    public static double getTrueOr(boolean[][] mask1, boolean[][] mask2) {
        int count= 0;
        for (int i= 0;i<mask1.length;i++) {
            for (int j= 0;j<mask1[0].length;j++) {
                if (mask1[i][j]==true||mask2[i][j]==true) count++;
            }   
        }
        return count;
    }
    
    public static double getTrueAndOr(boolean[][] mask1Or, boolean[][] mask2Or, boolean[] mask3And) {
        int count= 0;
        for (int i= 0;i<mask1Or.length;i++) {
            for (int j= 0;j<mask1Or[0].length;j++) {
                if (mask1Or[i][j]==true||mask2Or[i][j]==true) {
                    if (mask3And[i]==true) count++;
                        }
            }   
        }
        return count;
    }
    
    public static void makeMasks55k(Alignment known, Alignment unimputed, double cov, double het) {
        int covCutoff= (int)(cov*(double) unimputed.getSiteCount()); //cutoff for number of sites present to be high cov
        int hetCutoff= (int)(het*(double) unimputed.getSiteCount());
        int inbredCutoff= (int)(.006*(double) unimputed.getSiteCount());
        knownHetMask= new boolean[unimputed.getSequenceCount()][unimputed.getSiteCount()];
        calledHomoMajMask= new boolean[unimputed.getSequenceCount()][unimputed.getSiteCount()];
        calledHomoMinMask= new boolean[unimputed.getSequenceCount()][unimputed.getSiteCount()];
        calledMissingMask= new boolean[unimputed.getSequenceCount()][unimputed.getSiteCount()];
        highCovTaxa= new boolean[unimputed.getSequenceCount()];
        highCovInbreds= new boolean[unimputed.getSequenceCount()];
        highCovHets= new boolean[unimputed.getSequenceCount()];
        outbred= new boolean[unimputed.getSequenceCount()];
        inbred= new boolean[unimputed.getSequenceCount()];
        for (int taxon = 0; taxon < unimputed.getSequenceCount(); taxon++) {
            if (unimputed.getTotalNotMissingForTaxon(taxon)>covCutoff) highCovTaxa[taxon]= true;
            if (unimputed.getHeterozygousCountForTaxon(taxon)>hetCutoff) outbred[taxon]= true;
            if (unimputed.getHeterozygousCountForTaxon(taxon)<inbredCutoff) inbred[taxon]= true;
            if (highCovTaxa[taxon]==true&&outbred[taxon]==true) highCovHets[taxon]= true;
            if (highCovTaxa[taxon]==true&&inbred[taxon]==true) highCovInbreds[taxon]= true;
        }
        
        matchTaxon= new int[unimputed.getSequenceCount()];//holds the index of corresponding taxon in known
        String[] knownNames= new String[known.getSequenceCount()];
        for (int taxon = 0; taxon < knownNames.length; taxon++) {
            knownNames[taxon]= known.getIdGroup().getIdentifier(taxon).getNameLevel(0);
        }
        
        for (int taxon = 0; taxon < matchTaxon.length; taxon++) {
            String unkName= unimputed.getIdGroup().getIdentifier(taxon).getNameLevel(0);
            matchTaxon[taxon]= Arrays.binarySearch(knownNames, unkName);
        }
        
        int[] knownPos= known.getPhysicalPositions();
        for (int site = 0; site < unimputed.getSiteCount(); site++) {
            int matchSite= known.getSiteOfPhysicalPosition(unimputed.getPositionInLocus(site), null);
            if (Arrays.binarySearch(knownPos, unimputed.getPositionInLocus(site))<0) continue;
            byte diploidMaj= AlignmentUtils.getDiploidValue(unimputed.getMajorAllele(site), unimputed.getMajorAllele(site));
            byte diploidMin= AlignmentUtils.getDiploidValue(unimputed.getMinorAllele(site), unimputed.getMinorAllele(site));
            
            for (int taxon = 0; taxon < unimputed.getSequenceCount(); taxon++) {
                if (known.isHeterozygous(matchTaxon[taxon], matchSite)==true) knownHetMask[taxon][site]= true;
                else if (known.getBase(matchTaxon[taxon], matchSite)==diploidMaj) calledHomoMajMask[taxon][site]= true;
                else if (known.getBase(matchTaxon[taxon], matchSite)==diploidMin) calledHomoMinMask[taxon][site]= true;
                else if (known.getBase(matchTaxon[taxon],matchSite)==diploidN) calledMissingMask[taxon][site]= true;
            }
        }
        System.out.println("Parameters:\nHigh Coverage Cutoff: "+cov+"\nOutbred Cutoff:"+het+"\nInbred Cutoff: .006"
                    +"\n\nNumber of siteTaxa considered:\n"+"\nknownHets: "+
                    getTrue(knownHetMask)+" (highHet/highHomo: "+
                    getTrue(knownHetMask, highCovHets)+"/"+
                    getTrue(knownHetMask, highCovInbreds)+")"+"\n"+
                    "calledHomoMaj"+getTrue(calledHomoMajMask)+" (highHet/highHomo: "+
                    getTrue(calledHomoMajMask, highCovHets)+"/"+
                    getTrue(calledHomoMajMask, highCovInbreds)+")"+"\n"+
                    "calledHomoMin"+getTrue(calledHomoMinMask)+" (highHet/highHomo: "+
                    getTrue(calledHomoMinMask, highCovHets)+"/"+
                    getTrue(calledHomoMinMask, highCovInbreds)+")"+"\n"+
                    "calledMissing"+getTrue(calledMissingMask)+"\n\n"+
                    "Number of sequences in each group:\nTotalNumSequences: "+known.getSequenceCount()+"\nhighCovTaxa"+getTrue(highCovTaxa)+"\n"+
                    "highCovInbred"+getTrue(highCovInbreds)+"\n"+
                    "highCovHets"+getTrue(highCovHets)+"\n"+
                    "inbred"+getTrue(inbred)+"\n"+
                    "outbred"+getTrue(outbred)+"\n");
    }
    
    public static void makeMasks(Alignment known, int sampleIntensity, double cov, double het) { //het, homo, missing, highCov
        int covCutoff= (int)(cov*(double) known.getSiteCount()); //cutoff for number of sites present to be high cov
        int hetCutoff= (int)(het*(double) known.getSiteCount());
        int inbredCutoff= (int)(.006*(double) known.getSiteCount());
        knownHetMask= new boolean[known.getSequenceCount()][known.getSiteCount()];
        calledHomoMajMask= new boolean[known.getSequenceCount()][known.getSiteCount()];
        calledHomoMinMask= new boolean[known.getSequenceCount()][known.getSiteCount()];
        calledMissingMask= new boolean[known.getSequenceCount()][known.getSiteCount()];
        highCovTaxa= new boolean[known.getSequenceCount()];
        highCovInbreds= new boolean[known.getSequenceCount()];
        highCovHets= new boolean[known.getSequenceCount()];
        outbred= new boolean[known.getSequenceCount()];
        inbred= new boolean[known.getSequenceCount()];
        
        for (int taxon= 0;taxon<known.getSequenceCount();taxon++) {
            if (known.getTotalNotMissingForTaxon(taxon)>covCutoff) highCovTaxa[taxon]= true;
            if (known.getHeterozygousCountForTaxon(taxon)>hetCutoff) outbred[taxon]= true;
            if (known.getHeterozygousCountForTaxon(taxon)<inbredCutoff) inbred[taxon]= true;
            if (highCovTaxa[taxon]==true&&outbred[taxon]==true) highCovHets[taxon]= true;
            if (highCovTaxa[taxon]==true&&inbred[taxon]==true) highCovInbreds[taxon]= true;
            for (int site= taxon;site<known.getSiteCount();site+= sampleIntensity) {
                if (known.getBase(taxon, site)==diploidN) 
                    calledMissingMask[taxon][site]= true;
                else if (known.isHeterozygous(taxon, site)) 
                    knownHetMask[taxon][site]= true;                
                else if (known.getBaseArray(taxon, site)[0]==known.getMajorAllele(site)||known.getBaseArray(taxon, site)[1]==known.getMajorAllele(site)) 
                    calledHomoMajMask[taxon][site]= true;
                else if (known.getBaseArray(taxon, site)[0]==known.getMinorAllele(site)||known.getBaseArray(taxon, site)[1]==known.getMinorAllele(site))
                    calledHomoMinMask[taxon][site]= true;                
            }
        }
//        //get sites for high coverage hets in HW proportions
//        IdGroup highCovHetTaxa= IdGroupUtils.idGroupSubset(known.getIdGroup(), highCovHets);
//        Alignment highCovAlign= FilterAlignment.getInstance(known, highCovHetTaxa);
//        for (int site= 0;site<known.getSiteCount();site++) {
//            double p= highCovAlign.getMajorAlleleFrequency(site);
//            double q= highCovAlign.getMinorAlleleFrequency(site);
//            double obsHetFreq= highCovAlign.getHeterozygousCount(site)/highCovAlign.getSiteCount();
//            double expHetFreq= 2*p*q;
////            if (obsHetFreq>(expHetFreq-(q*hwWiggle))&&obsHetFreq<(expHetFreq+(q*hwWiggle))) HWSites[site]= true;
//        }
        //system.out to debug
        System.out.println("Parameters:\nHigh Coverage Cutoff: "+cov+"\nOutbred Cutoff:"+het+"\nInbred Cutoff: .006"
                    +"\nSample Intensity:"+sampleIntensity+"\n\nNumber of siteTaxa considered:\n"+"\nknownHets: "+
                    getTrue(knownHetMask)+" (highHet/highHomo: "+
                    getTrue(knownHetMask, highCovHets)+"/"+
                    getTrue(knownHetMask, highCovInbreds)+")"+"\n"+
                    "calledHomoMaj"+getTrue(calledHomoMajMask)+" (highHet/highHomo: "+
                    getTrue(calledHomoMajMask, highCovHets)+"/"+
                    getTrue(calledHomoMajMask, highCovInbreds)+")"+"\n"+
                    "calledHomoMin"+getTrue(calledHomoMinMask)+" (highHet/highHomo: "+
                    getTrue(calledHomoMinMask, highCovHets)+"/"+
                    getTrue(calledHomoMinMask, highCovInbreds)+")"+"\n"+
                    "calledMissing"+getTrue(calledMissingMask)+"\n\n"+
                    "Number of sequences in each group:\nTotalNumSequences: "+known.getSequenceCount()+"\nhighCovTaxa"+getTrue(highCovTaxa)+"\n"+
                    "highCovInbred"+getTrue(highCovInbreds)+"\n"+
                    "highCovHets"+getTrue(highCovHets)+"\n"+
                    "inbred"+getTrue(inbred)+"\n"+
                    "outbred"+getTrue(outbred)+"\n");
    }
    public static void QuickAccuracy(Alignment imputed, int taxon, int site) {
        if (knownBase!=diploidN&&impBase!=diploidN) {
            perSiteTaxon[knownIndex][imputed.getSequenceCount()]+=1.0;
            perSiteTaxon[knownSiteCount][taxon]+=1.0;
            if (AlignmentUtils.isEqual(knownBase, impBase)==false) {
//                System.out.println(NucleotideAlignmentConstants.getNucleotideIUPAC(knownBase)+"/"+NucleotideAlignmentConstants.getNucleotideIUPAC(impBase)+"site: "+knownIndex+"/"+site+"taxon: "+taxon);
                if (knownHet==true||imputed.isHeterozygous(taxon, site)==true) {
                    if (impArray[0]==knownArray[0]||impArray[1]==knownArray[0]||impArray[0]==knownArray[1]||
                        impArray[1]==knownArray[1]) perSiteTaxon[knownIndex][taxon]+=0;
                }
                else perSiteTaxon[knownIndex][taxon]+=1.0;
            }
                    
        }
    }
    
    public static void CalculateAccuracy(Alignment imputed, int taxon, int site, double[]type, double[]typeSites) {
        //calculates accuracy for all
                typeSites[0]++;
                if (calledHomoMajMask[taxon][site]==true) {
                    typeSites[1]++;
                    typeSites[2]++;
                    if (impBase==diploidN) type[5]++;
                    else if (knownBase==impBase) {
                        type[0]++;
                    }
                    else if (imputed.isHeterozygous(taxon, site)==true) {
                        if (impArray[0]==knownMaj||impArray[1]==knownMaj) type[2]++; 
                        else {
                            type[4]++;
                            type[14]++;
                        }
                    }
                    else {
                        type[3]++;
                        type[14]++;
                    }
                }
                else if (calledHomoMinMask[taxon][site]==true) {
                    typeSites[1]++;
                    typeSites[3]++;
                    if (impBase==diploidN) type[5]++;
                    else if (knownBase==impBase) type[1]++;
                    else if (imputed.isHeterozygous(taxon, site)==true) {
                        if (impArray[0]==knownMin||impArray[1]==knownMin) type[2]++;
                        else {
                            type[4]++;
                            type[14]++;
                        }
                    }
                    else {
                        type[3]++;
                        type[14]++;
                    }
                    
                }
                else if (knownHetMask[taxon][site]==true) {
                    typeSites[4]++;
                    if (impBase==diploidN) type[10]++;
                    else if (impBase==knownBase) type[6]++;
                    else if (imputed.isHeterozygous(taxon, site)==true) {
                        type[9]++;
                        type[14]++;
                    }
                    else if (impArray[0]==knownMaj||impArray[1]==knownMin||impArray[1]==knownMaj||impArray[0]==knownMin) {
                        type[7]++;
                        type[14]+= .5;
                    }
                    else {
                        type[8]++;
                        type[14]++;
                    }
                }
                else if (calledMissingMask[taxon][site]==true) {
                    typeSites[5]++;
                    if (impBase==diploidN) type[13]++;
                    else if (imputed.isHeterozygous(taxon, site)) type[12]++;
                    else type[11]++;
                }
    }
        
    public static void RunTest(Alignment known, Alignment imputed, Alignment unimputed, boolean knownTest, int sampleIntensity, double cov, double het, double hwWiggle, String outFileName) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        //for each double array, index 0:homoMajCorrect, 1:homoMinCorrect, 2:homoHetOneCorrect, 3:homoIncorrectHomo, 4:homoIncorrectHet
        //5:homoMissing, 6:hetCorrect, 7:hetOneCorrect, 8:hetIncorrectHomo, 9:hetIncorrectHet, 10:hetMissing, 11:missingImputedHomo, 
        //12:missingImputedHet, 13:missingMissing, 14:OverallError
        double[] all= new double[15];
        double[] allSites= new double[6]; //0 all, 1 homo, 2 homoMaj, 3 homoMin, 4 het, 5 missing
        double[] allInbred= new double[15];
        double[] allInbredSites= new double[6];
        double[] allOutbred= new double[15];
        double[] allOutbredSites= new double[6];
        double[] highCovInbred= new double[15];
        double[] highCovInbredSites= new double[6];
        double[] highCovOutbred= new double[15];
        double[] highCovOutbredSites= new double[6];
        knownSiteCount= known.getSiteCount();
        
        if (knownTest==true) {
            System.out.println("knownSitesMasked: "+known.getSiteCount()+"\nimputedTaxa: "+imputed.getSequenceCount());
            makeMasks55k(known, unimputed, cov, het);
            int[] knownPos= known.getPhysicalPositions();
            perSiteTaxon= new double[known.getSiteCount()+1][imputed.getSequenceCount()+1];
            for (int site = 0; site < imputed.getSiteCount(); site++) {
                knownIndex= Arrays.binarySearch(knownPos, imputed.getPositionInLocus(site));
                if (knownIndex<0) continue;
//                System.out.println(site+"/"+knownIndex);//debug
                int matchSite= known.getSiteOfPhysicalPosition(unimputed.getPositionInLocus(site), null);
                knownMaj= known.getMajorAllele(matchSite);
                knownMin= known.getMinorAllele(matchSite);
                for (int taxon = 0; taxon < unimputed.getSequenceCount(); taxon++) {
                    impArray= imputed.getBaseArray(taxon, site);
                    impBase= imputed.getBase(taxon, site);
                    knownBase= known.getBase(matchTaxon[taxon], matchSite);
                    knownArray= known.getBaseArray(matchTaxon[taxon], matchSite);
                    knownHet= known.isHeterozygous(matchTaxon[taxon], matchSite);
                    QuickAccuracy(imputed, taxon, site);
                    CalculateAccuracy(imputed,taxon,site,all,allSites);//calculates accuracy for all
                    if (inbred[taxon]==true) {//calculates accuracy for only those coded as inbred (not 1-het)
                        CalculateAccuracy(imputed,taxon,site,allInbred,allInbredSites);
                        if (highCovTaxa[taxon]==true) {//accuracy for inbred sites with high coverage
                            CalculateAccuracy(imputed,taxon,site,highCovInbred,highCovInbredSites);
                        }                    
                    }
                    if (outbred[taxon]==true) {//calculates accuracy for landraces only (heterozygosity above specified cutoff)
                        CalculateAccuracy(imputed,taxon,site,allOutbred,allOutbredSites);
                        if (highCovTaxa[taxon]==true) {//accuracy for landraces with high coverage
                            CalculateAccuracy(imputed,taxon,site,highCovOutbred,highCovOutbredSites);
                        }
                    }
                }
            }
        }
        else {
            makeMasks(known, sampleIntensity, cov, het);
            
            perSiteTaxon= new double[known.getSiteCount()+1][imputed.getSequenceCount()+1];
            for (int taxon= 0;taxon<known.getSequenceCount();taxon++) {            
                for (int site= taxon;site<known.getSiteCount();site+= sampleIntensity) {
                    knownIndex= site;
                    knownBase= known.getBase(taxon, site);
                    impBase= imputed.getBase(taxon, site);
                    knownMaj= known.getMajorAllele(site);
                    knownMin= known.getMinorAllele(site);
                    impArray= imputed.getBaseArray(taxon, site);
                    QuickAccuracy(imputed,taxon, site);
                    CalculateAccuracy(imputed,taxon,site,all,allSites);//calculates accuracy for all
                    if (inbred[taxon]==true) {//calculates accuracy for only those coded as inbred (not 1-het)
                        CalculateAccuracy(imputed,taxon,site,allInbred,allInbredSites);
                        if (highCovTaxa[taxon]==true) {//accuracy for inbred sites with high coverage
                            CalculateAccuracy(imputed,taxon,site,highCovInbred,highCovInbredSites);
                        }                    
                    }
                    if (outbred[taxon]==true) {//calculates accuracy for landraces only (heterozygosity above specified cutoff)
                        CalculateAccuracy(imputed,taxon,site,allOutbred,allOutbredSites);
                        if (highCovTaxa[taxon]==true) {//accuracy for landraces with high coverage
                            CalculateAccuracy(imputed,taxon,site,highCovOutbred,highCovOutbredSites);
                        }
                    }                
                }
            }
        }
        try{
            DecimalFormat df= new DecimalFormat("0.####");
            String outputFileName= dir+outFileName;
            DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
            
            outStream.writeBytes("\t\t\t\t\t\t\t\t\t\tcalledHomozygote\t\t\t\tcalledHeterozygote\t\t\t\t\tcalledMissing\n\tSitesConsidered\tSitesCalledHomo\tSiteCalledHomoMaj\tSitesCalledHomoMin\tSitesCalledHet\tSitesCalledMissing"
                    + "\thomoMajCorrect\thomoMinCorrect\thomoHetOneCorrect\thomoIncorrectHomo\thomoIncorrectHet"
                    + "\thomoMissing\thetCorrect\thetOneCorrect\thetIncorrectHomo\thetIncorrectHet\thetMissing\tmissingImputedHomo"
                    + "\tmissingImputedHet\tmissingMissing\tOverallError(halfHetsHalfCorrect)");
            outStream.writeBytes("\nallSites\t"+allSites[0]+"\t"+allSites[1]+"\t"+allSites[2]+"\t"+allSites[3]+"\t"+allSites[4]+"\t"+allSites[5]);
            outStream.writeBytes("\t"+all[0]/allSites[2]+"\t"+all[1]/allSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+all[i]/allSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+all[i]/allSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+all[i]/allSites[5]);}
            outStream.writeBytes("\t"+all[14]/(allSites[1]+allSites[4]));
            outStream.writeBytes("\nallInbred\t"+allInbredSites[0]+"\t"+allInbredSites[1]+"\t"+allInbredSites[2]+"\t"+allInbredSites[3]+"\t"+allInbredSites[4]+"\t"+allInbredSites[5]);
            outStream.writeBytes("\t"+allInbred[0]/allInbredSites[2]+"\t"+allInbred[1]/allInbredSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+allInbred[i]/allInbredSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+allInbred[i]/allInbredSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+allInbred[i]/allInbredSites[5]);}
            outStream.writeBytes("\t"+allInbred[14]/(allInbredSites[1]+allInbredSites[4]));
            outStream.writeBytes("\nhighCovInbred\t"+highCovInbredSites[0]+"\t"+highCovInbredSites[1]+"\t"+highCovInbredSites[2]+"\t"+highCovInbredSites[3]+"\t"+highCovInbredSites[4]+"\t"+highCovInbredSites[5]);
            outStream.writeBytes("\t"+highCovInbred[0]/highCovInbredSites[2]+"\t"+highCovInbred[1]/highCovInbredSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+highCovInbred[i]/highCovInbredSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+highCovInbred[i]/highCovInbredSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+highCovInbred[i]/highCovInbredSites[5]);}
            outStream.writeBytes("\t"+highCovInbred[14]/(highCovInbredSites[1]+highCovInbredSites[4]));
            outStream.writeBytes("\nallOutbred\t"+allOutbredSites[0]+"\t"+allOutbredSites[1]+"\t"+allOutbredSites[2]+"\t"+allOutbredSites[3]+"\t"+allOutbredSites[4]+"\t"+allOutbredSites[5]);
            outStream.writeBytes("\t"+allOutbred[0]/allOutbredSites[2]+"\t"+allOutbred[1]/allOutbredSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+allOutbred[i]/allOutbredSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+allOutbred[i]/allOutbredSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+allOutbred[i]/allOutbredSites[5]);}
            outStream.writeBytes("\t"+allOutbred[14]/(allOutbredSites[1]+allOutbredSites[4]));
            outStream.writeBytes("\nhighCovOutbred\t"+highCovOutbredSites[0]+"\t"+highCovOutbredSites[1]+"\t"+highCovOutbredSites[2]+"\t"+highCovOutbredSites[3]+"\t"+highCovOutbredSites[4]+"\t"+highCovOutbredSites[5]);
            outStream.writeBytes("\t"+highCovOutbred[0]/highCovOutbredSites[2]+"\t"+highCovOutbred[1]/highCovOutbredSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+highCovOutbred[i]/highCovOutbredSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+highCovOutbred[i]/highCovOutbredSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+highCovOutbred[i]/highCovOutbredSites[5]);}
            outStream.writeBytes("\t"+highCovOutbred[14]/(highCovOutbredSites[1]+highCovOutbredSites[4]));
            outStream.writeBytes("\ntaxonName:");
            for (int i= 0;i<imputed.getSequenceCount();i++) {outStream.writeBytes("\t"+imputed.getTaxaName(i));}
            outStream.writeBytes("\ntaxonSitesCompared:");
            for (int i= 0;i<perSiteTaxon[0].length;i++) {outStream.writeBytes("\t"+perSiteTaxon[known.getSiteCount()][i]);}
            outStream.writeBytes("\ntaxonError:");
            for (int taxon= 0;taxon<imputed.getSequenceCount();taxon++) {
                double bad= 0;
                for (int site= 0;site<known.getSiteCount();site++) {
                    bad+=perSiteTaxon[site][taxon];
                }
                double err= perSiteTaxon[known.getSiteCount()][taxon]!=0?((bad)/perSiteTaxon[known.getSiteCount()][taxon]):-1;
                outStream.writeBytes("\t"+df.format(err));
            }
            outStream.writeBytes("\nsiteName:");
            for (int i= 0;i<known.getSiteCount();i++) {outStream.writeBytes("\t"+known.getSNPID(i));}
            outStream.writeBytes("\nsiteSitesCompared:");
            for (int i= 0;i<perSiteTaxon.length;i++) {outStream.writeBytes("\t"+perSiteTaxon[i][imputed.getSequenceCount()]);}
            outStream.writeBytes("\nsiteError:");
            for (int site= 0;site<known.getSiteCount();site++) {
                double bad= 0;
                for (int taxon= 0;taxon<imputed.getSequenceCount();taxon++) {
                    bad+=perSiteTaxon[site][taxon];
                }
                double err= perSiteTaxon[site][imputed.getSequenceCount()]!=0?((bad)/perSiteTaxon[site][imputed.getSequenceCount()]):-1;
                outStream.writeBytes("\t"+df.format(err));
            }
            //to output the results matrix
//            for (int site = 0; site < perSiteTaxon.length; site++) {
//                for (int taxon = 0; taxon < perSiteTaxon[site].length; taxon++) {
//                    outStream.writeBytes(perSiteTaxon[site][taxon]+"\t");
//                }
//                outStream.writeBytes("\n");
//            }
            outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
    }
    
    public static void main(String[] args) {
        TasselPrefs.putAlignmentRetainRareAlleles(false);
        
        //make a mask
//        dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
//        String inFile= "SEED_12S_GBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1RndSample1000";
//        MaskFileSample(inFile,true,300);
        
//        dir= "/Users/kelly/Documents/GBS/Imputation/SmallFiles/";
//        MaskFile55k("RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.2", true,
//                "RIMMA_282_SNP55K_AGPv2_20100513__S45391.chr10_matchTo_RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1", true);
        
        //run accuracy
//        dir= "/Users/kelly/Documents/GBS/Imputation/";
////        dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
//        String knownFileName= "SEED_12S_GBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1";
//        String imputedFileName= "SEED_12S_GBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_masked_defaultDonor";
////        ImputationAccuracy.makeMasks(ImportUtils.readFromHapmap(dir+knownFileName+".hmp.txt", true, null), 300, .6, .02, .2);
//        ImputationAccuracy.RunTest(ImportUtils.readFromHapmap(dir+knownFileName+".hmp.txt.gz", null), 
//                ImportUtils.readFromHapmap(dir+imputedFileName+".hmp.txt.gz", null), null, false, 300,.6,.01,.2,imputedFileName+"Accuracy.6.txt");
        
        //run accuracy for 55k
//        dir= "/Users/kelly/Documents/GBS/Imputation/SmallFiles/";
//        ImputationAccuracy.RunTest(ImportUtils.readFromHapmap(dir+"RIMMA_282_SNP55K_AGPv2_20100513__S45391.chr10_matchTo_RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1.hmp.txt.gz",null),
//                    ImportUtils.readFromHapmap(dir+"RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.2_masked55k_HomoSegBlock200HetWithExtras8k.minMtCnt15.mxInbErr.01.mxHybErr.003.c10.hmp.txt",null),
//                    ImportUtils.readFromHapmap(dir+"RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.2.hmp.txt.gz",null),true, -1, .6, .005, .2, "RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.2_masked55k_HomoSegBlock200HetWithExtras8k.minMtCnt15.mxInbErr.01.mxHybErr.003.c10.Accuracy55k.txt");
        
        //mask input depth file using depth
//        dir= "/Users/kellyadm/Desktop/Imputation2.7/";//laptop
//        dir= "/home/local/MAIZE/kls283/GBS/Imputation2.7/parts/";
//        String h5Depth= dir+"AllZeaGBS_v2.7_SeqToGenos_part14.hmp.h5";
//        int depth= 5;
//        int maskDenom= 17;
//        maskBaseFileByDepth(h5Depth, depth, maskDenom, true);
//        
////         //mask using depth. requires an hdf5 file with depth
        dir= "/home/kls283/Documents/Imputation/";//cbsugbs
//        String h5Depth= dir+"AllZeaGBS_v2.7wDepth.hmp.h5";
        dir= "/Users/kellyadm/Desktop/Imputation2.7/";//laptop
//        String h5Depth= dir+"AllZeaGBS_v2.7_SeqToGenos_part14.hmp.h5";
//        String fileToMask= dir+"AllZeaGBS_v2.7_SeqToGenos_part14m.hmp.h5";
////        String fileToMask= dir+"AllZeaGBSv27StrictSubsetByAmes(no EP or GEM).hmp.h5";
////        String fileToMask= dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span.hmp.h5";
//        int depth= 5;
//        int maskDenom= 3;
////        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, false);
////        fileToMask= dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span.hmp.h5";
////        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, false);
//        fileToMask= dir+"AllZeaGBS_v2.7wDepthm.hmp.h5";
//        depth= 5;
//        maskDenom= 11;
//        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, true);
//        
//        fileToMask= dir+"AllZeaGBS_v2.7wDepthw.hmp.h5";
//        depth= 5;
//        maskDenom= 17;
//        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, true);
        
//        depth= 5;
//        maskDenom= 17;
//        fileToMask= dir+"AllZeaGBSv27StrictSubsetByAmes(no EP or GEM).hmp.h5";
//        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, false);
//        fileToMask= dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span.hmp.h5";
//        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, false);
//        
//        depth= 5;
//        maskDenom= 11;
//        fileToMask= dir+"AllZeaGBSv27StrictSubsetByAmes(no EP or GEM).hmp.h5";
//        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, false);
//        fileToMask= dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span.hmp.h5";
//        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, false);
//        
//        depth= 5;
//        maskDenom= 11;
//        fileToMask= dir+"AllZeaGBSv27StrictSubsetByAmes380.hmp.h5";
//        maskFileByDepth(h5Depth, fileToMask, depth, maskDenom,true, false);

        
        //run accuracy on depth-masked imputed file. Only system out now.
        String dir= "//Users/kls283/Desktop/Imputation/";
////        String dir= "/home/local/MAIZE/kls283/GBS/Imputation2.7/";
////        String key= dir+"AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span._maskKey_Depth5_Denom17.hmp.h5";
////        String imputed=  dir+"/results0918/AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span_LandraceDonor8k_depth5_denom17_imp.minCnt20.mxInbErr.01.mxHybErr.003.hmp.h5";
////        accuracyDepth(key, imputed);
////        imputed=  dir+"/results0918/AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span_StdDonor8k_depth5_denom17_imp.minCnt20.mxInbErr.01.mxHybErr.003.hmp.h5";
////        accuracyDepth(key, imputed);
////        imputed=  dir+"/results0918/AllZeaGBSv27StrictSubsetBy12S_RIMMA_Span_LandraceDonor4k_depth5_denom17_imp.minCnt20.mxInbErr.01.mxHybErr.003.hmp.h5";
////        accuracyDepth(key, imputed);
////        String key=  dir+"AllZeaGBSv27StrictSubsetByAmes408._maskKey_Depth5_Denom17.hmp.h5";
////        String imputed=  dir+"resultsBaseImputation/AllZeaGBSv27StrictSubsetByAmes408_LandraceDonor8k_depth5_denom17_imp.minCnt20.mxInbErr.01.mxHybErr.003.hmp.h5";
        String dataset= "NAM.rils.parents";//AmesTemperate 282All 12S_RIMMA_Span ExPVPIowa481
        String imputed=  dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy"+dataset+"NoIndelsMinTCov0.1MinSCov0.1Poly.2kDonorFRFocusToMiss_imp.minCnt20.mxInbErr.01.mxHybErr.003.hmp.h5";
        String key=  dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7StrictSubsetBy"+dataset+".hmp.h5";
        String mafFile= dir+"donors/AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7Match"+dataset+"_HaplotypeStd8kMAF.txt";
        double[] mafClass= new double[]{0,.02,.05,.1,.2,.3,.4,.5,1}; 
////        accuracyDepth(key, imputed);
////        imputed=  dir+"results0918/AllZeaGBSv27StrictSubsetByAmes408_StdDonor8k_depth5_denom17_imp.minCnt20.mxInbErr.01.mxHybErr.003.hmp.h5";
////        imputed= dir+"AllZeaGBSv27StrictSubsetByAmes408_LandraceDonor8kKellyFRFixed_depth5_denom17_imp.minCnt20.mxInbErr.01.mxHybErr.003.hmp.h5";
//        accuracyDepth(key, imputed, null, mafFile, mafClass);
        String peter= dir+"AllZeaGBS_v2.7_NAM.rils.parents.imputed.peters.algorithm.hmp.h5";
        String unimputed= dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetByNAM.rils.parentsNoIndelsMinTCov0.0MinSCov0.0Poly.hmp.h5";
//        accuracyDepth(key, peter,unimputed,null , mafClass);
        Alignment laura= ImportUtils.readFromHapmap(dir+"NAM.rils.laura.imp.hmp.txt", false, null);
        ExportUtils.writeToPlink(laura, dir+"NAM.rils.laura.imp.ped", '\t');
////        //for beagle
//        String dir= "//Users/kls283/Desktop/Imputation/";
//        String dataset= "SynB";//AmesTemperate 282All 12S_RIMMA_Span ExPVPIowa481 SynB
////        String chr= "8";
//////        String unimputed= dir+"beagle/AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy"+dataset+"chr"+chr+"NoIndels.vcf.gz";
//////        String imputed= dir+"beagle/AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy"+dataset+"chr"+chr+"NoIndelsBeagleBase.vcf.gz";
//////        String key=  dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7StrictSubsetBy"+dataset+"chr"+chr+".hmp.h5";
//        String unimputed= dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy"+dataset+"NoIndelsMinTCov0.1MinSCov0.1Poly.hmp.txt.gz";
//        String imputed= dir+"beagle/AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy"+dataset+"NoIndelsMinTCov0.1MinSCov0.1PolyBeagleBase.vcf.gz";
//        String key=  dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7StrictSubsetBy"+dataset+".hmp.h5";
//////        key= dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7StrictSubsetBy12S_RIMMA_Spanchr8.hmp.h5";
//        double[] MAFClass= new double[]{0,.02,.05,.1,.2,.3,.4,.5,1};//should always have a 0 (monomorphic) and 1 (>.5, shows that minor and major flipped)
//        accuracyDepth(key, imputed, unimputed, null, MAFClass);
        
//        //dummy imputation
//        dir= "//Users/kls283/Desktop/Imputation/";
//        String dataset= "SynB";//ExPVPIowa481 282All 12S_RIMMA_Span
//        String fileToImpute= dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy"+dataset+"NoIndelsMinTCov0.1MinSCov0.1Poly.hmp.txt.gz";
//        String keyFile= dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7StrictSubsetBy"+dataset+".hmp.h5";
//        imputeUsingPQ(fileToImpute);
//        accuracyDepth(keyFile,fileToImpute.substring(0, fileToImpute.indexOf(".hmp"))+"ImputedPQ.hmp.h5", null, null, null);
        
//        //joint accuracy
//        String dir= "//Users/kls283/Desktop/Imputation/";
//        String dataset= "SynB";//AmesTemperate 282All 12S_RIMMA_Span ExPVPIowa481
//        String imputed= dir+"AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy"+dataset+"NoIndelsMinTCov0.1MinSCov0.1Poly.4kDonorFRFocusToMiss_imp.minCnt20.mxInbErr.01.mxHybErr.003.hmp.h5";
//        String beagle= dir+"beagle/AllZeaGBS_v2.7wDepth_masked_Depth7_Denom7StrictSubsetBy"+dataset+"NoIndelsMinTCov0.1MinSCov0.1PolyBeagleBase.vcf.gz";
//        String key=  dir+"AllZeaGBS_v2.7wDepth_maskKey_Depth7_Denom7StrictSubsetBy"+dataset+".hmp.h5";
//        jointAccuracyDepth(key, imputed, beagle);
        
        }
    
}