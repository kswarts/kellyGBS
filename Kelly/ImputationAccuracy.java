package Kelly;

import java.io.*;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.popgen.MinorWindowViterbiImputation;

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
    public static byte diploidN= NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
    
    public static void maskFile(String inFile, boolean gz, int sampleIntensity) {
        Alignment a= ImportUtils.readFromHapmap(dir+inFile+(gz==true?".hmp.txt.gz":".hmp.txt"), null);
        MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
        for (int taxon= 0;taxon<a.getSequenceCount();taxon++) {            
            for (int site= taxon;site<a.getSiteCount();site+= sampleIntensity) {
                mna.setBase(taxon, site, diploidN);
            }
        }
        mna.clean();
        ExportUtils.writeToHapmap(mna, true, dir+inFile+"_masked.hmp.txt.gz", '\t', null);
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
        double[] allSites= new double[6]; //0 all, 1 homo, 2 homoMaj, 3 homoMin, 4 het, 5 missing
        double[] allInbred= new double[14];
        double[] allInbredSites= new double[6];
        double[] allOutbred= new double[14];
        double[] allOutbredSites= new double[6];
        double[] highCovInbred= new double[14];
        double[] highCovInbredSites= new double[6];
        double[] highCovOutbred= new double[14];
        double[] highCovOutbredSites= new double[6];
        
        ImputationAccuracy.makeMasks(known, sampleIntensity, cov, het, hwWiggle);
        for (int taxon= 0;taxon<known.getSequenceCount();taxon++) {            
            for (int site= taxon;site<known.getSiteCount();site+= sampleIntensity) {
                byte knownBase= known.getBase(taxon, site);
                byte impBase= imputed.getBase(taxon, site);
                byte[] knownArray= known.getBaseArray(taxon, site);
                byte knownMaj= known.getMajorAllele(site);
                byte knownMin= known.getMinorAllele(site);
                byte[] impArray= imputed.getBaseArray(taxon, site);
                //calculates accuracy for all
                allSites[0]++;
                if (calledHomoMajMask[taxon][site]==true) {
                    allSites[1]++;
                    allSites[2]++;
                    if (impBase==diploidN) all[5]++;
                    else if (knownBase==impBase) all[0]++;
                    else if (imputed.isHeterozygous(taxon, site)==true) {
                        if (impArray[0]==knownMaj||impArray[1]==knownMaj) all[2]++; 
                        else all[4]++;
                    }
                    else all[3]++;
                }
                else if (calledHomoMinMask[taxon][site]==true) {
                    allSites[1]++;
                    allSites[3]++;
                    if (impBase==diploidN) all[5]++;
                    else if (knownBase==impBase) all[1]++;
                    else if (imputed.isHeterozygous(taxon, site)==true) {
                        if (impArray[0]==knownMin||impArray[1]==knownMin) all[2]++;
                        else all[4]++;
                    }
                    else all[3]++;
                    
                }
                else if (knownHetMask[taxon][site]==true) {
                    allSites[4]++;
                    if (impBase==diploidN) all[10]++;
                    else if (impBase==knownBase) all[6]++;
                    else if (imputed.isHeterozygous(taxon, site)==true) all[9]++;
                    else if (impArray[0]==knownMaj||impArray[1]==knownMin||impArray[1]==knownMaj||impArray[0]==knownMin) all[7]++;
                    else all[8]++;
                }
                else if (calledMissingMask[taxon][site]==true) {
                    allSites[5]++;
                    if (impBase==diploidN) all[13]++;
                    else if (imputed.isHeterozygous(taxon, site)) all[12]++;
                    else all[11]++;
                }
                //calculates accuracy for only those coded as inbred (not 1-het)
                if (inbred[taxon]==true) {
                    allInbredSites[0]++;
                    if (calledHomoMajMask[taxon][site]==true) {
                        allInbredSites[1]++;
                        allInbredSites[2]++;
                        if (impBase==diploidN) allInbred[5]++;
                        else if (knownBase==impBase) allInbred[0]++;
                        else if (imputed.isHeterozygous(taxon, site)==true) {
                            if (impArray[0]==knownMaj||impArray[1]==knownMaj) allInbred[2]++;
                            else allInbred[4]++;
                        }
                        else allInbred[3]++;
                    }
                    else if (calledHomoMinMask[taxon][site]==true) {
                        allInbredSites[1]++;
                        allInbredSites[3]++;
                        if (impBase==diploidN) allInbred[5]++;
                        else if (knownBase==impBase) allInbred[1]++;
                        else if (imputed.isHeterozygous(taxon, site)==true) {
                            if (impArray[0]==knownMin||impArray[1]==knownMin) allInbred[2]++;
                            else allInbred[4]++;
                        }
                        else allInbred[3]++;

                    }
                    else if (knownHetMask[taxon][site]==true) {
                        allInbredSites[4]++;
                        if (impBase==diploidN) allInbred[10]++;
                        else if (impBase==knownBase) allInbred[6]++;
                        else if (imputed.isHeterozygous(taxon, site)==true) allInbred[9]++;
                        else if (impArray[0]==knownMaj||impArray[1]==knownMin||impArray[1]==knownMaj||impArray[0]==knownMin) allInbred[7]++;
                        else allInbred[8]++;
                    }
                    else if (calledMissingMask[taxon][site]==true) {
                        allInbredSites[5]++;
                        if (impBase==diploidN) allInbred[13]++;
                        else if (imputed.isHeterozygous(taxon, site)) allInbred[12]++;
                        else allInbred[11]++;
                    }
                    //accuracy for inbred sites with high coverage
                    if (highCovTaxa[taxon]==true) {
                        highCovInbredSites[0]++;
                        if (calledHomoMajMask[taxon][site]==true) {
                            highCovInbredSites[1]++;
                            highCovInbredSites[2]++;
                            if (impBase==diploidN) highCovInbred[5]++;
                            else if (knownBase==impBase) highCovInbred[0]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownMaj||impArray[1]==knownMaj) highCovInbred[2]++;
                                else highCovInbred[4]++;
                            }
                            else highCovInbred[3]++;
                        }
                        else if (calledHomoMinMask[taxon][site]==true) {
                            highCovInbredSites[1]++;
                            highCovInbredSites[3]++;
                            if (impBase==diploidN) highCovInbred[5]++;
                            else if (knownBase==impBase) highCovInbred[1]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownMin||impArray[1]==knownMin) highCovInbred[2]++; 
                                else highCovInbred[4]++;
                            }
                            else highCovInbred[3]++;                    
                        }
                        else if (knownHetMask[taxon][site]==true) {
                            highCovInbredSites[4]++;
                            if (impBase==diploidN) highCovInbred[10]++;
                            else if (impBase==knownBase) highCovInbred[6]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) highCovInbred[9]++;
                            else if (impArray[0]==knownMaj||impArray[1]==knownMin||impArray[1]==knownMaj||impArray[0]==knownMin) highCovInbred[7]++;
                            else highCovInbred[8]++;
                        }
                        else if (calledMissingMask[taxon][site]==true) {
                            highCovInbredSites[5]++;
                            if (impBase==diploidN) highCovInbred[13]++;
                            else if (imputed.isHeterozygous(taxon, site)) highCovInbred[12]++;
                            else highCovInbred[11]++;
                        }
                    }                    
                }
                //calculates accuracy for landraces only (heterozygosity above specified cutoff)
                if (outbred[taxon]==true) {
                    allOutbredSites[0]++;
                    if (calledHomoMajMask[taxon][site]==true) {
                        allOutbredSites[1]++;
                        allOutbredSites[2]++;
                        if (impBase==diploidN) allOutbred[5]++;
                        else if (knownBase==impBase) allOutbred[0]++;
                        else if (imputed.isHeterozygous(taxon, site)==true) {
                            if (impArray[0]==knownMaj||impArray[1]==knownMaj) allOutbred[2]++; 
                            else allOutbred[4]++;
                        }
                        else allOutbred[3]++;
                    }
                    else if (calledHomoMinMask[taxon][site]==true) {
                        allOutbredSites[1]++;
                        allOutbredSites[3]++;
                        if (impBase==diploidN) allOutbred[5]++;
                        else if (knownBase==impBase) allOutbred[1]++;
                        else if (imputed.isHeterozygous(taxon, site)==true) {
                            if (impArray[0]==knownMin||impArray[1]==knownMin) allOutbred[2]++; 
                            else allOutbred[4]++;
                        }
                        else allOutbred[3]++;

                    }
                    else if (knownHetMask[taxon][site]==true) {
                        allOutbredSites[4]++;
                        if (impBase==diploidN) allOutbred[10]++;
                        else if (impBase==knownBase) allOutbred[6]++;
                        else if (imputed.isHeterozygous(taxon, site)==true) allOutbred[9]++;
                        else if (impArray[0]==knownMaj||impArray[1]==knownMin||impArray[1]==knownMaj||impArray[0]==knownMin) allOutbred[7]++;
                        else allOutbred[8]++;
                    }
                    else if (calledMissingMask[taxon][site]==true) {
                        allOutbredSites[5]++;
                        if (impBase==diploidN) allOutbred[13]++;
                        else if (imputed.isHeterozygous(taxon, site)) allOutbred[12]++;
                        else allOutbred[11]++;
                    }
                    //accuracy for landraces with high coverage
                    if (highCovTaxa[taxon]==true) {
                        highCovOutbredSites[0]++;
                        if (calledHomoMajMask[taxon][site]==true) {
                            highCovOutbredSites[1]++;
                            highCovOutbredSites[2]++;
                            if (impBase==diploidN) highCovOutbred[5]++;
                            else if (knownBase==impBase) highCovOutbred[0]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownMaj||impArray[1]==knownMaj) highCovOutbred[2]++; 
                                else highCovOutbred[4]++;
                            }
                            else highCovOutbred[3]++;
                        }
                        else if (calledHomoMinMask[taxon][site]==true) {
                            highCovOutbredSites[1]++;
                            highCovOutbredSites[3]++;
                            if (impBase==diploidN) highCovOutbred[5]++;
                            else if (knownBase==impBase) highCovOutbred[1]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) {
                                if (impArray[0]==knownMin||impArray[1]==knownMin) highCovOutbred[2]++; 
                                else highCovOutbred[4]++;
                            }
                            else highCovOutbred[3]++;                    
                        }
                        else if (knownHetMask[taxon][site]==true) {
                            highCovOutbredSites[4]++;
                            if (impBase==diploidN) highCovOutbred[10]++;
                            else if (impBase==knownBase) highCovOutbred[6]++;
                            else if (imputed.isHeterozygous(taxon, site)==true) highCovOutbred[9]++;
                            else if (impArray[0]==knownMaj||impArray[1]==knownMin||impArray[1]==knownMaj||impArray[0]==knownMin) highCovOutbred[7]++;
                            else highCovOutbred[8]++;
                        }
                        else if (calledMissingMask[taxon][site]==true) {
                            highCovOutbredSites[5]++;
                            if (impBase==diploidN) highCovOutbred[13]++;
                            else if (imputed.isHeterozygous(taxon, site)) highCovOutbred[12]++;
                            else highCovOutbred[11]++;
                        }
                    }
                }                
            }
        }
        try{
            String outputFileName= dir+outFileName+"_accuracyTest.txt";
            DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outputFileName)));
            outStream.writeBytes("Parameters:\nHigh Coverage Cutoff: "+cov+"\nOutbred Cutoff:"+het+"\nInbred Cutoff: .006"
                    +"\nSample Intensity:"+sampleIntensity+"\n\nNumber of siteTaxa considered:\n"+"total number of sites: "+allSites[0]+"\nknownHets: "+
                    ImputationAccuracy.getTrue(knownHetMask)+" (highHet/highHomo: "+
                    ImputationAccuracy.getTrue(knownHetMask, highCovHets)+"/"+
                    ImputationAccuracy.getTrue(knownHetMask, highCovInbreds)+")"+"\n"+
                    "calledHomoMaj"+ImputationAccuracy.getTrue(calledHomoMajMask)+" (highHet/highHomo: "+
                    ImputationAccuracy.getTrue(calledHomoMajMask, highCovHets)+"/"+
                    ImputationAccuracy.getTrue(calledHomoMajMask, highCovInbreds)+")"+"\n"+
                    "calledHomoMin"+ImputationAccuracy.getTrue(calledHomoMinMask)+" (highHet/highHomo: "+
                    ImputationAccuracy.getTrue(calledHomoMinMask, highCovHets)+"/"+
                    ImputationAccuracy.getTrue(calledHomoMinMask, highCovInbreds)+")"+"\n"+
                    "calledMissing"+ImputationAccuracy.getTrue(calledMissingMask)+"\n\n"+
                    "Number of sequences in each group:\nTotalNumSequences: "+known.getSequenceCount()+"\nhighCovTaxa"+ImputationAccuracy.getTrue(highCovTaxa)+"\n"+
                    "highCovInbred"+ImputationAccuracy.getTrue(highCovInbreds)+"\n"+
                    "highCovHets"+ImputationAccuracy.getTrue(highCovHets)+"\n"+
                    "inbred"+ImputationAccuracy.getTrue(inbred)+"\n"+
                    "outbred"+ImputationAccuracy.getTrue(outbred)+"\n");
            outStream.writeBytes("\t\t\t\t\t\t\t\t\t\tcalledHomozygote\t\t\t\tcalledHeterozygote\t\t\t\t\t\tcalledMissing\n\tSitesConsidered\tSitesCalledHomo\tSiteCalledHomoMaj\tSitesCalledHomoMin\tSitesCalledHet\tSitesCalledMissing\thomoMajCorrect\thomoMinCorrect\thomoHetOneCorrect\thomoIncorrectHomo\thomoIncorrectHet"
                  + "\thomoMissing\thetCorrect\thetOneCorrect\thetIncorrectHomo\thetIncorrectHet\thetMissing\tmissingImputedHomo"
                  + "\tmissingImputedHet\tmissingMissing");
            outStream.writeBytes("\nallSites\t"+allSites[0]+"\t"+allSites[1]+"\t"+allSites[2]+"\t"+allSites[3]+"\t"+allSites[4]+"\t"+allSites[5]);
            outStream.writeBytes("\t"+all[0]/allSites[2]+"\t"+all[1]/allSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+all[i]/allSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+all[i]/allSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+all[i]/allSites[5]);}
            outStream.writeBytes("\nallInbred\t"+allInbredSites[0]+"\t"+allInbredSites[1]+"\t"+allInbredSites[2]+"\t"+allInbredSites[3]+"\t"+allInbredSites[4]+"\t"+allInbredSites[5]);
            outStream.writeBytes("\t"+allInbred[0]/allInbredSites[2]+"\t"+allInbred[1]/allInbredSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+allInbred[i]/allInbredSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+allInbred[i]/allInbredSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+allInbred[i]/allInbredSites[5]);}
            outStream.writeBytes("\nhighCovInbred\t"+highCovInbredSites[0]+"\t"+highCovInbredSites[1]+"\t"+highCovInbredSites[2]+"\t"+highCovInbredSites[3]+"\t"+highCovInbredSites[4]+"\t"+highCovInbredSites[5]);
            outStream.writeBytes("\t"+highCovInbred[0]/highCovInbredSites[2]+"\t"+highCovInbred[1]/highCovInbredSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+highCovInbred[i]/highCovInbredSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+highCovInbred[i]/highCovInbredSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+highCovInbred[i]/highCovInbredSites[5]);}
            outStream.writeBytes("\nallOutbred\t"+allOutbredSites[0]+"\t"+allOutbredSites[1]+"\t"+allOutbredSites[2]+"\t"+allOutbredSites[3]+"\t"+allOutbredSites[4]+"\t"+allOutbredSites[5]);
            outStream.writeBytes("\t"+allOutbred[0]/allOutbredSites[2]+"\t"+allOutbred[1]/allOutbredSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+allOutbred[i]/allOutbredSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+allOutbred[i]/allOutbredSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+allOutbred[i]/allOutbredSites[5]);}
            outStream.writeBytes("\nhighCovOutbred\t"+highCovOutbredSites[0]+"\t"+highCovOutbredSites[1]+"\t"+highCovOutbredSites[2]+"\t"+highCovOutbredSites[3]+"\t"+highCovOutbredSites[4]+"\t"+highCovOutbredSites[5]);
            outStream.writeBytes("\t"+highCovOutbred[0]/highCovOutbredSites[2]+"\t"+highCovOutbred[1]/highCovOutbredSites[3]);
            for (int i= 2;i<6;i++){outStream.writeBytes("\t"+highCovOutbred[i]/highCovOutbredSites[1]);}
            for (int i= 6;i<11;i++){outStream.writeBytes("\t"+highCovOutbred[i]/highCovOutbredSites[4]);}
            for (int i= 11;i<14;i++){outStream.writeBytes("\t"+highCovOutbred[i]/highCovOutbredSites[5]);}
            outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
    }
    
    public static void main(String[] args) {
        //run accuracy
//        dir= "/Users/kelly/Documents/GBS/GroupImputation/";
//        String knownFileName= "04_PivotMergedTaxaTBT.c10_s0_s24575";
//        String imputedFileName= "04_PivotMergedTaxaTBT.c10_s0_s24575_imputed";
////        ImputationAccuracy.makeMasks(ImportUtils.readFromHapmap(dir+knownFileName+".hmp.txt", true, null), 300, .6, .02, .2);
//        ImputationAccuracy.runTest(ImportUtils.readFromHapmap(dir+knownFileName+".hmp.txt", true, null), 
//                ImportUtils.readFromHapmap(dir+imputedFileName+".hmp.txt", true, null), 300, .6, .02, .2, imputedFileName);
        
        //make a mask
//        dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
//        String inFile= "SEED_12S_GBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1";
//        maskFile(inFile,true,300);
        
        //Ed's minorWindowViterbiImputation
        dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
        String donorFile= dir+"AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_HomoSegForTaxaGreaterThan.1Het_HaplotypeMerge_s+.hmp.txt.gz";
        String unImpTargetFile= dir+"SEED_12S_GBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_masked.hmp.txt.gz";
        String impTargetFile= dir+"SEED_12S_GBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_masked_imputedByNew.hmp.txt.gz";
        MinorWindowViterbiImputation e64NNI=new MinorWindowViterbiImputation(donorFile, unImpTargetFile, impTargetFile, 20, 50, 100, 0.01, true, false);
    }
}
