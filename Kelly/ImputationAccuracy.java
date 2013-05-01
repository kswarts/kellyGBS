package Kelly;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author kelly, alberto
 */
public class ImputationAccuracy {
    public static String dir;
    public static boolean[][] knownHetMask;
    public static boolean[][] calledHomoMajMask;
    public static boolean[][] calledHomoMinMask;
    public static boolean[][] calledMissingMask;
    public static boolean[] highCovTaxa;
    public static boolean[] highCovHets;
    public static boolean[] highCovInbreds;
    public static boolean[] HWSites;
    public static boolean[] outbred;
    public static boolean[] inbred;
    public static double[][] knownHetResults; //0 if mis/not imputed, 1 if one allele imputed correctly, -1 if both alleles imputed correctly
    public static double[][] calledHomoResults; //0 if mis/not imputed, 1 if allele imputed correctly, -1 if imputed as heterozygote
    public static double[][] calledMissingResults; //0 if not imputed, 1 if imputed as homozygote, 2 if imputed as heterozygote
    public static byte diploidN= NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");    
    
    public static int getTrue(boolean[] mask) {
        int count= 0;
        for (int i= 0;i<mask.length;i++) {
            if (mask[i]==true) count++;
        }
        return count;
    }  
    
    public static int getTrue(boolean[][] mask) {
        int count= 0;
        for (int i= 0;i<mask.length;i++) {
            for (int j= 0;j<mask[0].length;j++) {
                if (mask[i][j]==true) count++;
            }   
        }
        return count;
    }
    
    public static int getTrue(boolean[][] mask1, boolean[] mask2) {
        int count= 0;
        for (int i= 0;i<mask1.length;i++) {
            for (int j= 0;j<mask1[0].length;j++) {
                if (mask1[i][j]==true&&mask2[i]==true) count++;
            }   
        }
        return count;
    }
    
    public static void makeMasks(Alignment known, int sampleIntensity, double cov, double het, double hwWiggle) { //het, homo, missing, highCov, HW sites (wiggle refers to the leeway given as a proportion of the minor allele frequency)
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
        HWSites= new boolean[known.getSiteCount()];        
        
        for (int taxon= 0;taxon<known.getSequenceCount();taxon++) {
            if (known.getTotalNotMissingForTaxon(taxon)>covCutoff) highCovTaxa[taxon]= true;
            if (known.getHeterozygousCountForTaxon(taxon)>hetCutoff) outbred[taxon]= true;
            if (known.getHeterozygousCountForTaxon(taxon)<inbredCutoff) inbred[taxon]= true;
            if (highCovTaxa[taxon]==true&&outbred[taxon]==true) highCovHets[taxon]= true;
            if (highCovTaxa[taxon]==true&&inbred[taxon]==true) highCovInbreds[taxon]= true;
            for (int site= taxon;site<known.getSiteCount();site+= sampleIntensity) {
                if (known.getBase(taxon, site)==known.getMajorAllele(site)) 
                    calledHomoMajMask[taxon][site]= true;
                else if (known.getBase(taxon, site)==known.getMinorAllele(site))
                    calledHomoMinMask[taxon][site]= true;
                else if (known.getBase(taxon, site)==diploidN) 
                    calledMissingMask[taxon][site]= true;
                else if (known.isHeterozygous(taxon, site)) 
                    knownHetMask[taxon][site]= true;                
            }
        }
        //get sites for high coverage hets in HW proportions
        IdGroup highCovHetTaxa= IdGroupUtils.idGroupSubset(known.getIdGroup(), highCovHets);
        Alignment highCovAlign= FilterAlignment.getInstance(known, highCovHetTaxa);
        for (int site= 0;site<known.getSiteCount();site++) {
            double p= highCovAlign.getMajorAlleleFrequency(site);
            double q= highCovAlign.getMinorAlleleFrequency(site);
            double obsHetFreq= highCovAlign.getHeterozygousCount(site)/highCovAlign.getSiteCount();
            double expHetFreq= 2*p*q;
            if (obsHetFreq>(expHetFreq-(q*hwWiggle))&&obsHetFreq<(expHetFreq+(q*hwWiggle))) HWSites[site]= true;
        }
        //system.out to debug
        System.out.println("Number of siteTaxa considered:\n"+"knownHets: "+
                ImputationAccuracy.getTrue(knownHetMask)+" (highHet/highHomo: "+
                ImputationAccuracy.getTrue(knownHetMask, highCovHets)+"/"+
                ImputationAccuracy.getTrue(knownHetMask, highCovInbreds)+")"+"\n"+
                "calledHomoMaj"+ImputationAccuracy.getTrue(calledHomoMajMask)+" (highHet/highHomo: "+
                ImputationAccuracy.getTrue(calledHomoMajMask, highCovHets)+"/"+
                ImputationAccuracy.getTrue(calledHomoMajMask, highCovInbreds)+")"+"\n"+
                "calledHomoMin"+ImputationAccuracy.getTrue(calledHomoMinMask)+" (highHet/highHomo: "+
                ImputationAccuracy.getTrue(calledHomoMinMask, highCovHets)+"/"+
                ImputationAccuracy.getTrue(calledHomoMinMask, highCovInbreds)+")"+"\n"+
                "calledMissing"+ImputationAccuracy.getTrue(calledMissingMask)+"\n"+
                "Number of sequences in each group:\n"+"highCovTaxa"+ImputationAccuracy.getTrue(highCovTaxa)+"\n"+
                "highCovInbred"+ImputationAccuracy.getTrue(highCovInbreds)+"\n"+
                "highCovHets"+ImputationAccuracy.getTrue(highCovHets)+"\n"+
                "inbred"+ImputationAccuracy.getTrue(inbred)+"\n"+
                "outbred"+ImputationAccuracy.getTrue(outbred)+"\n"+
                "HWSites"+ImputationAccuracy.getTrue(HWSites));
    }
        
    public static void runTest(Alignment known, Alignment imputed, int sampleIntensity, double cov, double het, double hwWiggle, String outFileName) {
        //for each double array, index 0:homoMajCorrect, 1:homoMinCorrect, 2:homoHetOneCorrect, 3:homoIncorrectHomo, 4:homoIncorrectHet
        //5:homoMissing, 6:hetCorrect, 7:hetOneCorrect, 8:hetIncorrectHomo, 9:hetIncorrectHet, 10:hetMissing, 11:missingImputedHomo, 
        //12:missingImputedHet, 13:missingMissing
        double[] all= new double[14];
        double allSites= 0;
        double[] allInbred= new double[14];
        double allInbredSites= 0;
        double[] allOutbred= new double[14];
        double allOutbredSites= 0;
        double[] highCovInbred= new double[14];
        double highCovInbredSites= 0;
        double[] highCovOutbred= new double[14];
        double highCovOutbredSites= 0;
        
        ImputationAccuracy.makeMasks(known, sampleIntensity, cov, het, hwWiggle);
        for (int taxon= 0;taxon<known.getSequenceCount();taxon++) {            
            for (int site= taxon;site<known.getSiteCount();site+= sampleIntensity) {
                byte knownBase= known.getBase(taxon, site);
                byte impBase= imputed.getBase(taxon, site);
                byte[] knownArray= imputed.getBaseArray(taxon, site);
                byte[] impArray= imputed.getBaseArray(taxon, site);
                //calculates accuracy for all
                allSites++;
                if (calledHomoMajMask[taxon][site]==true) {
                    if (impBase==diploidN) all[5]++;
                    else if (knownBase==impBase) all[0]++;
                    else if (imputed.isHeterozygous(taxon, site)==true) {
                        if (impArray[0]==knownBase||impArray[1]==knownBase) all[2]++; 
                        else all[4]++;
                    }
                    else all[3]++;
                }
                else if (calledHomoMinMask[taxon][site]==true) {
                    if (impBase==diploidN) all[5]++;
                    else if (knownBase==impBase) all[1]++;
                    else if (imputed.isHeterozygous(taxon, site)==true) {
                        if (impArray[0]==knownBase||impArray[1]==knownBase) all[2]++; 
                        else all[4]++;
                    }
                    else all[3]++;
                    
                }
                else if (knownHetMask[taxon][site]==true) {
                    if (impBase==diploidN) all[10]++;
                    else if (impBase==knownBase) all[6]++;
                    else if (imputed.isHeterozygous(taxon, site)==true) all[9]++;
                    else if (impBase==knownArray[0]||impBase==knownArray[1]) all[7]++;
                    else all[8]++;
                }
                else if (calledMissingMask[taxon][site]==true) {
                    if (impBase==diploidN) all[13]++;
                    else if (imputed.isHeterozygous(taxon, site)) all[12]++;
                    else all[11]++;
                }
                //calculates accuracy for only those coded as inbred (not 1-het)
                if (inbred[taxon]==true) {
                    //accuracy for inbred sites with high coverage
                    if (highCovTaxa[taxon]==true) {
                        highCovInbredSites++;
                        if (calledHomoMajMask[taxon][site]==true) {
                            if (impBase==diploidN) highCovInbred[5]++;
                            else if (knownBase==impBase) highCovInbred[0]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownBase||impArray[1]==knownBase) highCovInbred[2]++; 
                                else highCovInbred[4]++;
                            }
                            else highCovInbred[3]++;
                        }
                        else if (calledHomoMinMask[taxon][site]==true) {
                            if (impBase==diploidN) highCovInbred[5]++;
                            else if (knownBase==impBase) highCovInbred[1]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownBase||impArray[1]==knownBase) highCovInbred[2]++; 
                                else highCovInbred[4]++;
                            }
                            else highCovInbred[3]++;                    
                        }
                        else if (knownHetMask[taxon][site]==true) {
                            if (impBase==diploidN) highCovInbred[10]++;
                            else if (impBase==knownBase) highCovInbred[6]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) highCovInbred[9]++;
                            else if (impBase==knownArray[0]||impBase==knownArray[1]) highCovInbred[7]++;
                            else highCovInbred[8]++;
                        }
                        else if (calledMissingMask[taxon][site]==true) {
                            if (impBase==diploidN) highCovInbred[13]++;
                            else if (imputed.isHeterozygous(taxon, site)) highCovInbred[12]++;
                            else highCovInbred[11]++;
                        }
                    }
                    //code accuracy for all sites coded as inbred
                    else {
                        allInbredSites++;
                        if (calledHomoMajMask[taxon][site]==true) {
                            if (impBase==diploidN) allInbred[5]++;
                            else if (knownBase==impBase) allInbred[0]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownBase||impArray[1]==knownBase) allInbred[2]++; 
                                else allInbred[4]++;
                            }
                            else allInbred[3]++;
                        }
                        else if (calledHomoMinMask[taxon][site]==true) {
                            if (impBase==diploidN) allInbred[5]++;
                            else if (knownBase==impBase) allInbred[1]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownBase||impArray[1]==knownBase) allInbred[2]++; 
                                else allInbred[4]++;
                            }
                            else allInbred[3]++;                    
                        }
                        else if (knownHetMask[taxon][site]==true) {
                            if (impBase==diploidN) allInbred[10]++;
                            else if (impBase==knownBase) allInbred[6]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) allInbred[9]++;
                            else if (impBase==knownArray[0]||impBase==knownArray[1]) allInbred[7]++;
                            else allInbred[8]++;
                        }
                        else if (calledMissingMask[taxon][site]==true) {
                            if (impBase==diploidN) allInbred[13]++;
                            else if (imputed.isHeterozygous(taxon, site)) allInbred[12]++;
                            else allInbred[11]++;
                        }
                    }
                }
                //calculates accuracy for landraces only (heterozygosity above specified cutoff)
                if (outbred[taxon]==true) {
                    //accuracy for landraces with high coverage
                    if (highCovTaxa[taxon]==true) {
                        highCovOutbredSites++;
                        if (calledHomoMajMask[taxon][site]==true) {
                            if (impBase==diploidN) highCovOutbred[5]++;
                            else if (knownBase==impBase) highCovOutbred[0]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownBase||impArray[1]==knownBase) highCovOutbred[2]++; 
                                else highCovOutbred[4]++;
                            }
                            else highCovOutbred[3]++;
                        }
                        else if (calledHomoMinMask[taxon][site]==true) {
                            if (impBase==diploidN) highCovOutbred[5]++;
                            else if (knownBase==impBase) highCovOutbred[1]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownBase||impArray[1]==knownBase) highCovOutbred[2]++; 
                                else highCovOutbred[4]++;
                            }
                            else highCovOutbred[3]++;                    
                        }
                        else if (knownHetMask[taxon][site]==true) {
                            if (impBase==diploidN) highCovOutbred[10]++;
                            else if (impBase==knownBase) highCovOutbred[6]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) highCovOutbred[9]++;
                            else if (impBase==knownArray[0]||impBase==knownArray[1]) highCovOutbred[7]++;
                            else highCovOutbred[8]++;
                        }
                        else if (calledMissingMask[taxon][site]==true) {
                            if (impBase==diploidN) highCovOutbred[13]++;
                            else if (imputed.isHeterozygous(taxon, site)) highCovOutbred[12]++;
                            else highCovOutbred[11]++;
                        }
                    }
                    else {
                        allOutbredSites++;
                        if (calledHomoMajMask[taxon][site]==true) {
                            if (impBase==diploidN) allOutbred[5]++;
                            else if (knownBase==impBase) allOutbred[0]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownBase||impArray[1]==knownBase) allOutbred[2]++; 
                                else allOutbred[4]++;
                            }
                            else allOutbred[3]++;
                        }
                        else if (calledHomoMinMask[taxon][site]==true) {
                            if (impBase==diploidN) allOutbred[5]++;
                            else if (knownBase==impBase) allOutbred[1]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownBase||impArray[1]==knownBase) allOutbred[2]++; 
                                else allOutbred[4]++;
                            }
                            else allOutbred[3]++;                    
                        }
                        else if (knownHetMask[taxon][site]==true) {
                            if (impBase==diploidN) allOutbred[10]++;
                            else if (impBase==knownBase) allOutbred[6]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) allOutbred[9]++;
                            else if (impBase==knownArray[0]||impBase==knownArray[1]) allOutbred[7]++;
                            else allOutbred[8]++;
                        }
                        else if (calledMissingMask[taxon][site]==true) {
                            if (impBase==diploidN) allOutbred[13]++;
                            else if (imputed.isHeterozygous(taxon, site)) allOutbred[12]++;
                            else allOutbred[11]++;
                        }
                    }
                }                
            }
        }
        try{
            String outFile= dir+outFileName+"_accuracyTest.txt";
          DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile)));
          outStream.writeBytes("\thomoMajCorrect\thomoMinCorrect\thomoHetOneCorrect\thomoIncorrectHomo\thomoIncorrectHet"
                  + "\thomoMissing\thetCorrect\thetOneCorrect\thetIncorrectHomo\thetIncorrectHet\thetMissing\tmissingImputedHomo"
                  + "\tmissingImputedHet\tmissingMissing");
          outStream.writeBytes("\nallSites");
          for (int i= 0;i<all.length;i++){outStream.writeBytes("\t"+all[i]/allSites);}
          outStream.writeBytes("\nallInbred");
          for (int i= 0;i<allInbred.length;i++){outStream.writeBytes("\t"+allInbred[i]/allInbredSites);}
          outStream.writeBytes("\nhighCovInbred");
          for (int i= 0;i<highCovInbred.length;i++){outStream.writeBytes("\t"+highCovInbred[i]/highCovInbredSites);}
          outStream.writeBytes("\nallOutbred");
          for (int i= 0;i<allInbred.length;i++){outStream.writeBytes("\t"+allOutbred[i]/allOutbredSites);}
          outStream.writeBytes("\nhighCovOutbred");
          for (int i= 0;i<highCovOutbred.length;i++){outStream.writeBytes("\t"+highCovOutbred[i]/highCovOutbredSites);}
          outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
    }
    
    public static void main(String[] args) {
        dir= "/Users/kelly/Documents/GBS/GroupImputation/";
        String knownFileName= "04_PivotMergedTaxaTBT.c10_s0_s24575";
        String imputedFileName= "04_PivotMergedTaxaTBT.c10_s0_s24575_imputed";
//        ImputationAccuracy.makeMasks(ImportUtils.readFromHapmap(dir+knownFileName+".hmp.txt", true, null), 300, .6, .02, .2);
        ImputationAccuracy.runTest(ImportUtils.readFromHapmap(dir+knownFileName+".hmp.txt", true, null), 
                ImportUtils.readFromHapmap(dir+imputedFileName+".hmp.txt", true, null), 300, .6, .02, .2, dir+imputedFileName);
        
    }
}
