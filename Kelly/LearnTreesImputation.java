package Kelly;


import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.text.AttributedCharacterIterator.Attribute;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import net.maizegenetics.gbs.pipeline.MergeIdenticalGametes;
import net.maizegenetics.pal.alignment.*;
import net.maizegenetics.pal.distance.IBSDistanceMatrix;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.time.StopWatch;
import sun.security.jca.GetInstance.Instance;
import weka.core.*;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 * Learning algorithms for imputation
 * @author kelly
 */
public class LearnTreesImputation {
    public static byte haploidN= 0xf;
    public static byte diploidN= NucleotideAlignmentConstants.getNucleotideDiploidByte("NN");
    public static int startIndex;
    public static int endIndex;    
    public StopWatch time= new StopWatch();
    
//    public static String[] GetGenotypeAtSite(Alignment a, int site) {
//        String[] alleles= new String[a.getSequenceCount()];
//        byte base;
//        for (int taxon=0; taxon<a.getSequenceCount();taxon++) {
//            base= a.getBase(taxon,site);
//            alleles[taxon]= NucleotideAlignmentConstants.NUCLEOTIDE_DIPLOID_HASH.get(base);
//        }
//        return alleles;
//    }
    
        public static String[] GetIUPACAtSite(Alignment a, int site) {
        String[] alleles= new String[a.getSequenceCount()];
        for (int taxon=0; taxon<a.getSequenceCount();taxon++) {
            alleles[taxon]= NucleotideAlignmentConstants.getNucleotideIUPAC(a.getBase(taxon, site));
        }
        return alleles;
    }
    
    public static byte[][] GetAllelesTwoSites(Alignment a, int site1, int site2) { //returns a byte 2-D array of alleles at the input sites
        byte[][] allelesInByte= new byte[2][a.getSequenceCount()];
        for (int i=0; i<a.getSequenceCount();i++) {
            allelesInByte[0][i]= a.getBase(i, site1);
            allelesInByte[1][i]= a.getBase(i, site2);
        }
        
        return allelesInByte;
    }
    
    public static int GetMissing(Alignment a, int site) {
        int missing= 0;
        for (int taxon= 0;taxon<a.getSequenceCount();taxon++) {
            if (a.getBase(taxon,site)==diploidN) missing++; 
        }
        return missing;
    }
    
    public static int NumPolymorphicSites(Alignment a) {
        int poly= 0;
        for (int s= 0; s<a.getSiteCount(); s++) {
            if (a.isPolymorphic(s) == true) poly++;
        }
        return poly;
    }
    
    public static boolean[] MaskAlleles(Alignment a, int site, byte haploidAllele) {//mask a site based on an allele in byte (haploid value), haploidAllele value in constructor is true
        boolean[] mask= new boolean[a.getSequenceCount()];
        byte baseArray[];
        for (int i=0; i<mask.length; i++){
            baseArray= a.getBaseArray(i,site);
            if (baseArray[0] == haploidAllele||baseArray[1] == haploidAllele) mask[i]= true;
        }
        
        return mask;
    }
    
    public static Alignment SelectForTaxaPresent(Alignment a, int site) { //returns a new alignment filtered first for sample, then for polymorphic sites
        //filter for taxa
        IdGroup IDs= a.getIdGroup();
        boolean[] includeTaxon= new boolean[a.getSequenceCount()];
        Arrays.fill(includeTaxon,true);
        for (int taxa= 0; taxa<a.getSequenceCount(); taxa++) {
            if (a.getBase(taxa,site)!=diploidN) includeTaxon[taxa]= false;
        }
        IdGroup removeIDs= IdGroupUtils.idGroupSubset(IDs, includeTaxon);
        Alignment align= FilterAlignment.getInstanceRemoveIDs(a, removeIDs);
        //filter for sites that are polymorphic
        ArrayList<Integer> subSite= new ArrayList<Integer>();
        for (int s= 0; s<align.getSiteCount(); s++) {
            if (align.isPolymorphic(s) == true) subSite.add(s);
        }
        int[] keepSite= ArrayUtils.toPrimitive(subSite.toArray(new Integer[LearnTreesImputation.NumPolymorphicSites(align)]));
        Alignment finalAlign= FilterAlignment.getInstance(align, keepSite);
        return finalAlign;
    }
    
    public static Alignment SelectForTaxaAbsent(Alignment a, int site) {
        //filter for taxa
        IdGroup IDs= a.getIdGroup();
        boolean[] includeTaxon= new boolean[a.getSequenceCount()];
        Arrays.fill(includeTaxon,true);
        for (int taxa= 0; taxa<a.getSequenceCount(); taxa++) {
            if (a.getBase(taxa,site)==diploidN) includeTaxon[taxa]= false;
        }
        IdGroup removeIDs= IdGroupUtils.idGroupSubset(IDs, includeTaxon);
        Alignment align= FilterAlignment.getInstanceRemoveIDs(a, removeIDs);
        
        return align;
    }
    
    public static Alignment SubsetAlignmentByAllele(Alignment a, int site, byte haploidAllele) {//make a subset alignment that only contains sequences where a site contains a haploid value specified
        IdGroup IDs= a.getIdGroup();
        boolean[] include= new boolean[a.getSequenceCount()];
        Arrays.fill(include,true);
        byte[] baseArray;
        for (int taxon= 0; taxon<a.getSequenceCount(); taxon++) {
            baseArray= a.getBaseArray(taxon, site);
            if (baseArray[0] == haploidAllele || baseArray[1] == haploidAllele) include[taxon]= false;
        }
        IdGroup removeIDs= IdGroupUtils.idGroupSubset(IDs, include);
        Alignment align= FilterAlignment.getInstanceRemoveIDs(a, removeIDs);    
        
        return align;
    }
    
    public static double GetEntropy(Alignment a, int siteIndex) {
        double Pmaj= a.getMajorAlleleFrequency(siteIndex);
        double Pmin= a.getMinorAlleleFrequency(siteIndex);
        double Pmiss= 1-((double) a.getTotalGametesNotMissing(siteIndex)/((double) a.getSequenceCount()*2));
        double Smaj= Pmaj==0.0?0.0:-Math.log10(Pmaj);
        double Smin= Pmin==0.0?0.0:-Math.log10(Pmin);
        double Smiss= Pmiss==0.0?0.0:-Math.log10(Pmiss);
        double H= (Pmaj*Smaj)+(Pmin*Smin)*(Pmiss*Smiss);
        
        return H;
    }
    
    public static double GetConditionalEntropy(Alignment a, int site2, boolean[] mask) { //mask must be same length as getSequenceCount(); site2 is the site being compared to the reference site; the mask is for the haploid allele value of the reference site that conditions the entropy
        double gametes= 0;
        for (int i=0;i<mask.length;i++){if (mask[i]==true) gametes+= 2.0;}
        double Pmaj= a.getMajorAlleleFrequencyForSubset(site2, mask);
        double Pmin= a.getMinorAlleleFrequencyForSubset(site2, mask);
        double Pmiss= 1-(((double) a.getTotalGametesNotMissingForSubset(site2, mask))/gametes);
        double Smaj= Pmaj==0.0?0.0:-Math.log10(Pmaj);
        double Smin= Pmin==0.0?0.0:-Math.log10(Pmin);
        double Smiss= Pmiss==0.0?0.0:-Math.log10(Pmiss);
        double H= (Pmaj*Smaj)+(Pmin*Smin)*(Pmiss*Smiss);
        
        return H;
    }
                
    
    public static double[] CalculateLDForSite(Alignment a, int siteIndex, int window) {
        Alignment align= FilterAlignment.getInstance(a, startIndex, endIndex);
        LinkageDisequilibrium theLD = new LinkageDisequilibrium(align, -1,LinkageDisequilibrium.testDesign.SiteByAll, siteIndex, null, false, -1, null);
        theLD.run();
        double[] LD= new double[align.getSiteCount()];
        for (int site2= 0;site2<align.getSiteCount();site2++) {
            LD[site2]= theLD.getRSqr(siteIndex, site2);            
        }
        return LD;
    }
    
    public static double CalculateGeneticDistanceForGroupBasedOnHaploidAllele(Alignment a, int siteIndex, byte haploidAllele, boolean forMissing) {
        Alignment align;
        if (forMissing==false) align= LearnTreesImputation.SubsetAlignmentByAllele(a, siteIndex, haploidAllele);
        else align= LearnTreesImputation.SelectForTaxaAbsent(a, siteIndex);
        IBSDistanceMatrix dm= new IBSDistanceMatrix(align);
        double avg= dm.meanDistance();
        
        return avg;
    }
    
    public static int FindBestSite(Alignment a, int siteIndex, int window, String out) { //calculates the information gain with respect to the target site for each other site, returns index of site with greatest gain        
        double maxInfo= 0;
        int bestSiteIndex= siteIndex;
        double Pmaj= a.getMajorAlleleFrequency(siteIndex);
        byte majAllele= a.getMajorAllele(siteIndex);
        double Pmin= a.getMinorAlleleFrequency(siteIndex);
        byte minAllele= a.getMinorAllele(siteIndex);
        double[] infoGain= new double[a.getSiteCount()];
        startIndex= (siteIndex-window<0)?0:siteIndex-window;
        endIndex= (siteIndex+window>a.getSiteCount())?a.getSiteCount()-1:siteIndex+window;
//        double[] LD= LearnTreesImputation.CalculateLDForSite(a,siteIndex,window);
//        System.out.println("index site: "+a.getSNPID(siteIndex)+"\n");
        
        for (int site2= startIndex;site2<endIndex;site2++) {
            infoGain[site2]= Pmaj*LearnTreesImputation.GetConditionalEntropy(a,site2,LearnTreesImputation.MaskAlleles(a,siteIndex,majAllele))
                    +Pmin*LearnTreesImputation.GetConditionalEntropy(a,site2,LearnTreesImputation.MaskAlleles(a,siteIndex,minAllele));
            if (infoGain[site2]>maxInfo) {
                bestSiteIndex= site2;
                maxInfo= infoGain[site2];
            }
//            System.out.println(infoGain[site2]);
        }
        try{
            String outFile= "/Users/kelly/Documents/GBS/TestImputation/FindBestSiteTest_"+out+"_bestSite"+bestSiteIndex+"new.txt";
          DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 655360));
          outStream.writeBytes("SiteIndex\tSNPID\tInfoGain\tSiteCov\tHetCount\tMAF\tPhysicalDist");
          
          for (int site= startIndex; site < endIndex; site++) {              
              outStream.writeBytes("\n"+site+"\t"+a.getSNPID(site)+"\t"+infoGain[site]+"\t"+(double)a.getTotalGametesNotMissing(site)/((double)a.getSequenceCount()*2)+"\t"+a.getHeterozygousCount(site) +"\t"+a.getMinorAlleleFrequency(site)+"\t"+Math.abs(a.getPositionInLocus(siteIndex)-a.getPositionInLocus(site)));//+"\t"+LD);              
           }
          outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
        
        return bestSiteIndex;
    }
    
//    public static int[][] GetTreeForSite(Alignment a, int refSite) {
//        int[][] tree;
//        
//        return tree;
//    }
    
    public static Alignment GetHaplotypes(Alignment a, int windowSize) {
        
    }
    
    public static Instances CreateWekaInstances(Alignment a, String instancesName) {
        ArrayList<Attribute> atts= new ArrayList<Attribute>(a.getSiteCount());
        String[] values= new String[3];
        values[2]= "N";
        for(int site= 0; site<a.getSiteCount();site++) {
            values[0]= a.getMajorAlleleAsString(site);
            values[1]= a.getMinorAlleleAsString(site);
            List<String> valueList= Arrays.asList(values);
            atts.add(site,new Attribute(a.getSNPID(site), valueList));
        }
        Instances input= new Instances(instancesName, atts, a.getSequenceCount());
        for (int taxon= 0;taxon<a.getSequenceCount();taxon++) {
            String[] vals= new String[a.getSiteCount()];
            for (int site= 0;site<a.getSiteCount();site++) {
                vals[site]= a.getBaseAsString(taxon, site);
            }
            input.add(taxon, new Instance(vals));
        }
        return input;
    }
    
    public static int[][] BuildTree(Alignment a, int siteIndex, int window) {
        int[][] tree= new int[a.getSequenceCount()][a.getSequenceCount()];
        startIndex= (siteIndex-window<0)?0:siteIndex-window;
        endIndex= (siteIndex+window>a.getSiteCount())?a.getSiteCount()-1:siteIndex+window;
        
        return tree;
    }
    
    public static void PruneTree() { //prune based on accuracy in validation set (maybe add complexity penalty?)
        
    }
    
    public static void main (String args[]) {
       //for get best site, doesn't actually generate a tree
//       String dir= "/Users/kelly/Documents/GBS/WGSHapmap/";
//       String inFile= "maizeHapMapV2_B73RefGenV2_201203028_chr8subset_129000000-135000000Polymorphic";
//       Alignment a= ImportUtils.readFromHapmap(dir+inFile+".hmp.txt", true, null);  
//       System.out.println("Samples in hapmap: "+a.getSequenceCount()+"\n"+"Sites in hapmap: "+a.getSiteCount());
//       int refSite= 50000;
//       int window= 100000;
//       String SNP= a.getSNPID(refSite);
//       Alignment noMissing= LearnTreesImputation.SelectForTaxaPresent(a, refSite);
//       String[] allelesPre= LearnTreesImputation.GetIUPACAtSite(a, refSite);       
//       System.out.println("Samples to train on: "+noMissing.getSequenceCount()+"\n"+"Sites remaining after filtering: "+noMissing.getSiteCount());
//       double H= LearnTreesImputation.GetEntropy(noMissing, refSite);
//       System.out.println("H for "+SNP+": "+H);
//       for (int i=0;i<allelesPre.length;i++){System.out.println(allelesPre[i]);}
//       int root= LearnTreesImputation.FindBestSite(noMissing, refSite, window, inFile+"_site"+refSite+"_window"+window);
//       
//       double[] GDRoot= new double[3];
//       GDRoot[0]= LearnTreesImputation.CalculateGeneticDistanceForGroupBasedOnHaploidAllele(a, root, a.getMajorAllele(root), false);
//       GDRoot[1]= LearnTreesImputation.CalculateGeneticDistanceForGroupBasedOnHaploidAllele(a, root, a.getMinorAllele(root), false);
//       GDRoot[2]= LearnTreesImputation.CalculateGeneticDistanceForGroupBasedOnHaploidAllele(a, root, a.getMajorAllele(root), true);
//       System.out.println("Best Site: "+noMissing.getSNPID(root)+" ("+root+")"+"\n"+"GD for major allele: "+GDRoot[0]+"\n"+"GD for minor allele: "+GDRoot[1]+"\n"+"GD for missing: "+GDRoot[2]);
        
   }
}
