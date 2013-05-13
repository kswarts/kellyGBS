/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import static Kelly.KellyUtils.dir;
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
import org.apache.commons.lang.ArrayUtils;

/**
 *
 * @author kls283
 */
public class AssignHaplotypes {
    public static String dir;
    public static byte diploidN= (byte) 0xff;
    
    public static void MatchSitesToAlignment(String inFileRef, boolean gzRef, String inFileMod, boolean gzMod) {
       String inFileRefName= (gzRef==true)?dir+inFileRef+".hmp.txt.gz":dir+inFileRef+".hmp.txt";
       String inFileModName= (gzMod==true)?dir+inFileMod+".hmp.txt.gz":dir+inFileMod+".hmp.txt";
       String outFileName= dir+inFileMod+"_sitesMatch"+inFileRef+".hmp.txt";
       Alignment ref= ImportUtils.readFromHapmap(inFileRefName, null);
       Alignment mod= ImportUtils.readFromHapmap(inFileModName, null);
       ArrayList<Integer> subSite= new ArrayList<Integer>();
       int[] refPos= ref.getPhysicalPositions();
       int currModPos= 0;
       //check to make sure that maj/min are the same between alignments for physical positions that match
       int sitesWithSamePos= 0;
       int disagree= 0;
       try{
           DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(dir+inFileMod+inFileRef+".SystemOutput.txt"), 655360));
           outStream.writeBytes("Maj/min allele of sites in modified/reference file that do not match site at corresponding physical position:");
           for (int site= 0;site<mod.getSiteCount();site++) {
               currModPos= mod.getPositionInLocus(site);
               int refIndex= Arrays.binarySearch(refPos, currModPos);
               if (refIndex>0) {
                   sitesWithSamePos++;
                   if ((ref.getMajorAllele(refIndex)==mod.getMajorAllele(site)&&ref.getMinorAllele(refIndex)==mod.getMinorAllele(site))||
                           (ref.getMajorAllele(refIndex)==mod.getMinorAllele(site)&&ref.getMinorAllele(refIndex)==mod.getMajorAllele(site))) subSite.add(site);
                   else {
                       outStream.writeBytes("Physical position: "+mod.getPositionInLocus(site)+"\tSiteIndex: "+site+"/"+refIndex+"\tMod/Ref Maj: ("+mod.getMajorAlleleAsString(site)+"/"+ref.getMajorAlleleAsString(refIndex)+")"+"\tMod/Ref Min: ("+mod.getMinorAlleleAsString(site)+"/"+ref.getMinorAlleleAsString(refIndex)+")");
                       disagree++;
                   }
               }
           }
           outStream.writeBytes(disagree+" out of "+sitesWithSamePos+" sites with same physical position do not share the same maj/min allele");
           outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
       subSite.trimToSize();
       int[] keepSite= ArrayUtils.toPrimitive(subSite.toArray(new Integer[subSite.size()]));
       Alignment sub= FilterAlignment.getInstance(mod, keepSite);
       ExportUtils.writeToHapmap(sub, true, outFileName, '\t', null);
   }
    
    public static void AssignInbredIdentityToHaps(String hapFile, boolean gzHap, String inbredRefFile, boolean gzInbred, int siteBlock, double mismatchTol) {
       String hapFileName= (gzHap==true)?dir+hapFile+".hmp.txt.gz":dir+hapFile+".hmp.txt";
       String inbredFileName= (gzInbred==true)?dir+inbredRefFile+".hmp.txt.gz":dir+inbredRefFile+".hmp.txt";
       String outFileStringName= dir+hapFile+"MatchedForIdentityTo"+inbredFileName+"_String.txt";
       String outFileIndexName= dir+hapFile+"MatchedForIdentityTo"+inbredFileName+"_InbredIndexNum.txt";
       Alignment hap= ImportUtils.readFromHapmap(hapFileName, null);
       Alignment inbred= ImportUtils.readFromHapmap(inbredFileName, null);
       
       int[][] index= new int[hap.getSequenceCount()][(hap.getSiteCount()/siteBlock)+1]; //instatiate array to hold inbred assignments to -1
       for (int i=0;i<index.length;i++) {
           for (int j=0;j<index[0].length;j++) {
               index[i][j]= -1;
           }
       }
       //assigns the highest matching inbred taxon to haplotype chunk, assuming any meet a minimum cutoff. -1 is unassigned. assume complete inbreeding
       for (int block= 0; block<siteBlock; block++) {
           for (int haplo= 0; haplo<hap.getSequenceCount();haplo++) {
               double prevTaxonSim= 0.0;
               double currTaxonSimilarity= 0.0;
               for (int taxon= 0; taxon<inbred.getSequenceCount(); taxon++) {
                   int countSame= 0;
                   int comparedSites= 0;
                   for (int site= block+0; site<Math.min(block+siteBlock,hap.getSequenceCount()); site++) {
                       if (hap.getBase(taxon, site)!=diploidN && inbred.getBase(taxon, site)!=diploidN) {
                           comparedSites++;
                           if (hap.getBase(taxon, site)==inbred.getBase(taxon, site)) countSame++;
                       }
                   }
                   prevTaxonSim= currTaxonSimilarity;
                   currTaxonSimilarity=countSame/comparedSites;
                   if (prevTaxonSim<currTaxonSimilarity && currTaxonSimilarity>1-mismatchTol) index[haplo][block]= taxon;
               } 
           }
       }
       try{
            DataOutputStream outIndexStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFileIndexName), 655360));
            DataOutputStream outStringStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFileStringName), 655360));
            outIndexStream.writeBytes("haplotype file: "+hapFileName+"\ninbred file: "+inbredFileName);
            for (int i=0;i<index.length;i++) {
                outIndexStream.writeBytes("\n"+hap.getTaxaName(i)+"\t");
                outStringStream.writeBytes("\n"+hap.getTaxaName(i)+"\t");
                for (int j=0;j<index[0].length;j++) {
                    outIndexStream.writeBytes("\t"+index[i][j]);
                    outStringStream.writeBytes("\t"+inbred.getTaxaName(index[i][j]));
            }
        }
            outStringStream.close();
            outIndexStream.close();
        }
        
        catch(IOException e) {
           System.out.println(e);
        }
    }
    
    public static void main(String[] args) {
        //for matchSitesInAlignment
       dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
       String inMod= "04_PivotMergedTaxaTBT.c10_s0_s24575subset__minCov0.1";
       String inRef= "maizeHapMapV2_B73RefGenV2_201203028_chr10";
       MatchSitesToAlignment(inRef,false,inMod,false);
    }
}
