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
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
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
       String outFileName= dir+inFileRef+"_with"+inFileMod+"TaxaIncluded.hmp.txt.gz";
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
           outStream.writeBytes("RefFile: "+inFileRefName+"\nModFile: "+inFileModName+"\nMaj/min allele of sites in modified/reference file that do not match site at corresponding physical position:");
           for (int site= 0;site<mod.getSiteCount();site++) {
               currModPos= mod.getPositionInLocus(site);
               int refIndex= Arrays.binarySearch(refPos, currModPos);
               if (refIndex>-1) {
                   sitesWithSamePos++;
                   if ((ref.getMajorAllele(refIndex)==mod.getMajorAllele(site)&&ref.getMinorAllele(refIndex)==mod.getMinorAllele(site))||
                           (ref.getMajorAllele(refIndex)==mod.getMinorAllele(site)&&ref.getMinorAllele(refIndex)==mod.getMajorAllele(site))) subSite.add(site);
                   else {
                       outStream.writeBytes("\nPhysical position: "+mod.getPositionInLocus(site)+"\tSiteIndex: "+site+"/"+refIndex+"\tMod/Ref Maj: ("+mod.getMajorAlleleAsString(site)+"/"+ref.getMajorAlleleAsString(refIndex)+")"+"\tMod/Ref Min: ("+mod.getMinorAlleleAsString(site)+"/"+ref.getMinorAlleleAsString(refIndex)+")");
                       disagree++;
                   }
               }
           }
           outStream.writeBytes("\n"+disagree+" out of "+sitesWithSamePos+" sites with same physical position do not share the same maj/min allele");
           System.out.println(disagree+" out of "+sitesWithSamePos+" sites with same physical position do not share the same maj/min allele");
           outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
       subSite.trimToSize();
       int[] keepSite= ArrayUtils.toPrimitive(subSite.toArray(new Integer[subSite.size()]));
       Alignment sub= FilterAlignment.getInstance(mod, keepSite);
       
       //modify IdGroup for the mod file names to reflect origin
       SimpleIdGroup newNames=  SimpleIdGroup.getInstance(sub.getIdGroup());
       for (int name= 0;name<newNames.getIdCount();name++) {
           newNames.setIdentifier(name, Identifier.getMergedInstance(new Identifier("HappyMap"), newNames.getIdentifier(name)));
       }
       System.out.println("Generating new hapmap with "+sub.getSequenceCount()+" additional taxa with information at "+sub.getSiteCount()+" sites");
       int[] subPos= sub.getPhysicalPositions();
       int currRefPos= 0;
       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(ref, mod.getSequenceCount()+ref.getSequenceCount(), ref.getSiteCount());
       for (int taxon=0;taxon<sub.getSequenceCount();taxon++) {
           mna.addTaxon(newNames.getIdentifier(taxon));
           System.out.println("Adding taxon: "+newNames.getName(taxon));
           for (int site= 0;site<ref.getSiteCount();site++) {
               currRefPos= ref.getPositionInLocus(site);
               int subIndex= Arrays.binarySearch(subPos, currRefPos);
               if (subIndex>-1) mna.setBase(taxon+ref.getSequenceCount(), site, sub.getBase(taxon, subIndex));
           }
       }
       
       mna.clean();
       ExportUtils.writeToHapmap(mna, true, outFileName, '\t', null);
   }
    
    public static void AssignInbredIdentityToHaps(String hapFile, boolean gzHap, String inbredRefFile, boolean gzInbred, int siteBlock, double mismatchTol) {
       String hapFileName= (gzHap==true)?dir+hapFile+".hmp.txt.gz":dir+hapFile+".hmp.txt";
       String inbredFileName= (gzInbred==true)?dir+inbredRefFile+".hmp.txt.gz":dir+inbredRefFile+".hmp.txt";
       String outTaxaStringName= dir+hapFile+"MatchedForIdentityTo"+inbredFileName+"_String.txt";
       String outTaxaIndexName= dir+hapFile+"MatchedForIdentityTo"+inbredFileName+"_InbredIndexNum.txt";
       String outHaploIndexName= dir+hapFile+"MatchedForIdentityTo"+inbredFileName+"_HaploIndexNum.txt";
       Alignment hap= ImportUtils.readFromHapmap(hapFileName, null);
       Alignment inbred= ImportUtils.readFromHapmap(inbredFileName, null);
       
       int[][] inbredIndex= new int[hap.getSequenceCount()][(hap.getSiteCount()/siteBlock)+1]; //instatiate array to hold inbred assignments to -1
       for (int i=0;i<inbredIndex.length;i++) {
           for (int j=0;j<inbredIndex[0].length;j++) {
               inbredIndex[i][j]= -1;
           }
       }
       int[][] hapIndex= new int[inbred.getSequenceCount()][(inbred.getSiteCount()/siteBlock)+1]; //instantiate array to hold haplotype assignments to inbred
       for (int i=0;i<hapIndex.length;i++) {
           for (int j=0;j<hapIndex[0].length;j++) {
               hapIndex[i][j]= -1;
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
                       if (hap.getBase(taxon, site)!=diploidN && inbred.getBase(taxon, site)!=diploidN && 
                               hap.isHeterozygous(taxon, site)==false && inbred.isHeterozygous(taxon, site)==false) {
                           comparedSites++;
                           if (hap.getBase(taxon, site)==inbred.getBase(taxon, site)) countSame++;
                       }
                   }
                   prevTaxonSim= currTaxonSimilarity;
                   currTaxonSimilarity=countSame/comparedSites;
                   if (prevTaxonSim<currTaxonSimilarity && currTaxonSimilarity>1-mismatchTol) {
                       inbredIndex[haplo][block]= taxon;
                       hapIndex[taxon][block]= haplo;
                   }
               } 
           }
       }
       try{
            DataOutputStream outTaxaIndexStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outTaxaIndexName), 655360));
            DataOutputStream outTaxaStringStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outTaxaStringName), 655360));
            DataOutputStream outHaploIndexStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outTaxaStringName), 655360));
            outTaxaIndexStream.writeBytes("haplotype file: "+hapFileName+"\ninbred file: "+inbredFileName);
            outTaxaStringStream.writeBytes("haplotype file: "+hapFileName+"\ninbred file: "+inbredFileName);
            outHaploIndexStream.writeBytes("haplotype file: "+hapFileName+"\ninbred file: "+inbredFileName);
            for (int i=0;i<inbredIndex.length;i++) {
                outTaxaIndexStream.writeBytes("\n"+hap.getTaxaName(i)+"\t");
                outTaxaStringStream.writeBytes("\n"+hap.getTaxaName(i)+"\t");
                for (int j=0;j<inbredIndex[0].length;j++) {
                    outTaxaIndexStream.writeBytes("\t"+inbredIndex[i][j]);
                    outTaxaStringStream.writeBytes("\t"+inbred.getTaxaName(inbredIndex[i][j]));
            }
        }
            for (int i=0;i<hapIndex.length;i++) {
                outHaploIndexStream.writeBytes("\n"+inbred.getTaxaName(i)+"\t");
                for (int j=0;j<inbredIndex[0].length;j++) {
                    outHaploIndexStream.writeBytes("\t"+hapIndex[i][j]);
            }
        }
            outTaxaStringStream.close();
            outTaxaIndexStream.close();
            outHaploIndexStream.close();
        }
        
        catch(IOException e) {
           System.out.println(e);
        }
    }
    
    public static void main(String[] args) {
        //for matchSitesInAlignment
       dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
       String inRef= "AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1";
       String inMod= "maizeHapMapV2_B73RefGenV2_201203028_chr10";
       MatchSitesToAlignment(inRef,false,inMod,false);
    }
}
