package Kelly;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */



import java.io.*;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Scanner;
import java.util.Set;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.gbs.util.BaseEncoder;
import net.maizegenetics.pal.alignment.Alignment;
import net.maizegenetics.pal.alignment.AlignmentUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.FilterAlignment;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignment;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pal.alignment.MutableNucleotideDepthAlignment;
import net.maizegenetics.pal.alignment.NucleotideAlignmentConstants;
import net.maizegenetics.pal.ids.IdGroup;
import net.maizegenetics.pal.ids.IdGroupUtils;
import net.maizegenetics.pal.ids.Identifier;
import net.maizegenetics.pal.ids.SimpleIdGroup;
import net.maizegenetics.pal.popgen.LinkageDisequilibrium;
import net.maizegenetics.prefs.TasselPrefs;
import net.maizegenetics.util.ExceptionUtils;
import net.maizegenetics.util.Utils;
import org.apache.commons.lang.ArrayUtils;
import org.apache.commons.lang.time.FastDateFormat;

/**
 *
 * @author kelly
 */
public class KellyUtils {
    public static String dir;
    public static Date theDate = new Date();
    public static FastDateFormat fdf = FastDateFormat.getInstance("yyyy-MM-dd:HH-mm-ss");
    public static String theDateString = fdf.format(theDate);
    public static byte diploidN= (byte) 0xff;

    
   public static void TagsByTaxaToFastq(String inFile, String outFile) {
      TagsByTaxaByte theTags= new TagsByTaxaByte(inFile, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type

      String tagSequence;
      int tagLength;
      int count;
      try{
          DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outFile), 655360));

          for (int tag= 0; tag < theTags.getTagCount(); tag++) {
              tagLength= theTags.getTagLength(tag);
              count= theTags.getReadCount(tag);
              outStream.writeBytes("@length="+tagLength+"count="+count+"\n");   //Length & count header
              tagSequence = BaseEncoder.getSequenceFromLong(theTags.getTag(tag));
              tagSequence = tagSequence.substring(0, tagLength);  //Remove any poly-A padding
              outStream.writeBytes(tagSequence+"\n+\n");    //Sequence and "+" symbol
              for (int i= 0; i < tagLength; i++) { outStream.writeBytes("f"); }           //Bogus quality string
                  outStream.writeBytes("\n");
           }
       }

      catch(IOException e) {
           System.out.println(e);
       }
   }

   //only implemented for tbt.Byte files just now. really only need to use if need to retain read counts for TagTaxon
   public static void FilterMergedTBT(int minSharedTaxa, String inTBTFileName, String outTBTFileName) {
       TagsByTaxaByte theTBT= new TagsByTaxaByte(inTBTFileName, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       TagCountMutable outputTags = new TagCountMutable(theTBT.getTag(0).length, theTBT.getTagCount()); //new tagsMutable the size of theTBT
       String[] theTaxa=theTBT.getTaxaNames();

       for (int tag= 0; tag < theTBT.getTagCount(); tag++) {//adds the tags that appear in greater than or equal to the taxon specified in the constructor
           if (theTBT.getNumberOfTaxaWithTag(tag) > minSharedTaxa) {
               outputTags.addReadCount(theTBT.getTag(tag), theTBT.getTagLength(tag), 1);
           }
       }
       outputTags.shrinkToCurrentRows();

       TagsByTaxaByte theTBTOut = new TagsByTaxaByte(theTaxa, outputTags);//new blank TBT.Byte file with the taxon from theTBT and the filtered tags
       for (int tagFill=0; tagFill < theTBTOut.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
           for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
               theTBTOut.addReadsToTagTaxon(tagFill, taxonFill, theTBT.getReadCountForTagTaxon(theTBT.getTagIndex(theTBTOut.getTag(tagFill)), taxonFill));
           }
       }
       theTBTOut.writeDistFile(new File(outTBTFileName), FilePacking.Byte, 1);
   }
   
//   public static void DissectTOPM(String inTOPMFileName) {
//       byte[] byteTest = new byte[1];
//       byteTest[0] = Byte.MIN_VALUE;
//       System.out.println(byteTest[0]);
//       File outTextFileName = new File(inTOPMFileName + ".txt");
//       TagsOnPhysicalMap theTOPM = new TagsOnPhysicalMap(inTOPMFileName, true);
//       int[] dist = theTOPM.mappingDistribution();
//
//       System.out.println("\nTotal tags:\t" + theTOPM.tagNum
//               + "\nMapped tags:\t" + theTOPM.mappedTags()              
//               + "\nTag Distribution:\n");
//       for (int dis = 0; dis < dist.length; dis++) {
//           System.out.println("Tags mapping to " + dis + " locations:\t" + dist[dis]);
//       }
//       theTOPM.writeTextFile(outTextFileName);
//
//
////       int noMap= 0;
////       int oneMap= 0;
////       int twoMap= 0;
////       int greaterTwoMap= 0;
////       a
////       for (int site= 0; site < theTOPM.tagNum; site++) {
////           theTOPM.getPositionArray(site);
////       }
//   }
   /**for biologically characterizing a TOPM and associated tagsByTaxaFile. Positions determined based on AGP_v1, centromeres taken from Wolfgruber et al 2012 and pericentromeric/telomeric extrapolated from Gore et al 2009**/ 
   public static void CharTOPMWithTBT(String inTOPMFileName, String inTBTFileName, int divisor) {
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads in the TOPM, true if binary
       TagsByTaxaByteHDF5TagGroups theTBT= new TagsByTaxaByteHDF5TagGroups(inTBTFileName);
       File outputFile= new File(dir+"/CharTOPMWithTBTLog_unitsOf"+divisor+"TagsFromTBT"+theDateString);
       
       //order TBT by read count
       int[][] orderByRead= new int[theTBT.getTagCount()][2];//a 2D int array to hold site indices from the topm and the corresponding read counts
       for (int sortTBT= 0; sortTBT < theTBT.getTagCount(); sortTBT++) {//fill orderedMasterTags with indices (index 1) and read counts (index 0)
           orderByRead[sortTBT][0]= theTBT.getReadCount(sortTBT);
           orderByRead[sortTBT][1]= sortTBT;
       }

       Arrays.sort(orderByRead, new java.util.Comparator<int[]>() {//sort 2d array by read count, descending
           @Override
           public int compare(int[] o1, int[] o2) {
               if(o1[0]>o2[0]) return -1; // -1 for descending order
               else if(o1[0]<o2[0]) return 1; // 1 for descending order
               else return 0;
           }
       });
       
       int[][] stats= new int[(theTBT.getTagCount()/divisor)+1][18];//index 0 (which division), 1 (alignCount), 2(multiMap), 3(telomeric), 4(pericentromeric), 5(centromeric, 6-17(chromosome count)
       double[][] taxaStats= new double[theTBT.getTagCount()/divisor+1][2];
       int[][] chrPos= {//where chr1 is index 0, chr2 is index 1 and so on 
               {133300000,89300000,94600000,104200000,101600000,49800000,55300000,45900000,68600000,59300000},//holds the starting centromeric position for each chromosome using the map position of functional centromeres from Wolfgruber at al Table 1 
               {133900000,91100000,95400000,105000000,108600000,50400000,55700000,48000000,69200000,60700000},//holds the ending centromeric position for each chromosome using the map position of functional centromeres from Wolfgruber at al Table 1. Use the greatest range for chr5
               {100000000,65000000,30000000,40000000,40000000,35000000,50000000,25000000,35000000,25000000},//holds starting "pericentromeric" position extrapolated from Gore et al 2009 based on where the recombination rate decreases
               {150000000,160000000,120000000,130000000,140000000,85000000,90000000,100000000,85000000,75000000}//holds ending "pericentromeric" position extrapolated from Gore et al 2009 based on where the recombination rate increases
           };
       
       int division= 0;
       int align= 0;
       int notAlign= 0;
       int taxaNumOfAlignedTags= 0;       
       int taxaNumOfUnalignedTags= 0;
       
       int realTBTIndex;
       long[] tag;
       int TOPMIndex;
       int chr;
       int minPos;
       
       for (int orderedTBTIndex= 0; orderedTBTIndex < theTBT.getTagCount(); orderedTBTIndex++) {
           realTBTIndex= orderByRead[orderedTBTIndex][1];
           tag= theTBT.getTag(realTBTIndex);
           TOPMIndex= theTOPM.getTagIndex(tag);
           chr=  theTOPM.getChromosome(TOPMIndex);
           minPos= theTOPM.getStartPosition(TOPMIndex);
                      
           if (theTOPM.getStartPosition(TOPMIndex) > -1) {
               align++;
               taxaNumOfAlignedTags+= theTBT.getNumberOfTaxaWithTag(realTBTIndex); //num taxon that share aligned site
               stats[division][chr+6]++;//plus 6 because 6 other variables held in array (for chromosomes 0-12)
               stats[division][1]++;//num of tags that align
               
               if (chr > 0 && chr < 11) {
                   if ((minPos < chrPos[1][chr-1]) && (minPos > chrPos[0][chr-1])) stats[division][5]++;                     
                   else if (minPos < chrPos[3][chr-1] && minPos > chrPos[2][chr-1]) stats[division][4]++;                     
                   else stats[division][3]++;
               }
           }
           else {
               taxaNumOfUnalignedTags+= theTBT.getNumberOfTaxaWithTag(realTBTIndex);//num taxon that share unaligned site
               notAlign++;
           }
           
           if (theTOPM.getPositionArray(TOPMIndex).length > 1) stats[division][2]++;
           
           if (((orderedTBTIndex)%divisor == 0 && orderedTBTIndex != 0) || theTBT.getTagCount()-1 == orderedTBTIndex) {
               taxaStats[division][0]= ((double)taxaNumOfAlignedTags/(double)(align*theTBT.getTaxaCount()))*(double)theTBT.getTaxaCount();//the average number of taxon that share an aligned site in the divisor set
               taxaStats[division][1]= ((double)taxaNumOfUnalignedTags/(double)(notAlign*theTBT.getTaxaCount()))*(double)theTBT.getTaxaCount();//the average number of taxon that share an unaligned site in the divisor set
               division++;
               stats[division-1][0]= division;
               taxaNumOfAlignedTags= 0;       
               taxaNumOfUnalignedTags= 0;
               align= 0;
               notAlign= 0;
               System.out.println("finish "+(division*divisor)+" tags in the TBT, ordered by read count");
               
           }
       }
       try {
                    outputFile.createNewFile();
                    outputFile.setWritable(true);
                    BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));

                    BW.write("Division\tAligned\tMultipleMapping\tNumAlignedTelomeric\tNumAlignedPericentromeric\tNumAlignedCentromeric");
                    for (int i= 0; i < 13; i++) {
                        BW.write("\tChr"+i);
                    }
                    BW.write("\tAvgNumTaxaWithAlignedTag\tAvgNumTaxaWithUnalignedTag\n");
                    
                    for (int div=0; div < stats.length; div++) {
                        for (int content= 0; content < 18; content++) {
                            BW.write(stats[div][content]+"\t");                            
                        }
                        BW.write(taxaStats[div][0]+"\t"+taxaStats[div][1]);
                        
                        BW.write("\n");
                    }
                    BW.close();
                }

                catch (IOException e) {
                    System.out.println(e);
                }
   }

    /**for biologically characterizing a TOPM and associated tagsByTaxaFile. Positions determined based on AGP_v1, centromeres taken from Wolfgruber et al 2012 and pericentromeric/telomeric extrapolated from Gore et al 2009**/ 
   public static void CharTOPMWithTagCount(String inTOPMFileName, String inTagCountFileName, int divisor) {
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads in the TOPM, true if binary
       TagCounts theTags= new TagCounts(inTagCountFileName, FilePacking.Byte);
       File outputFile= new File(dir+"/CharTOPMWithTagCountLog_unitsOf"+divisor+"TagsFromTBT"+theDateString);
       
       //order TBT by read count
       int[][] orderByRead= new int[theTags.getTagCount()][2];//a 2D int array to hold site indices from the topm and the corresponding read counts
       for (int sortTags= 0; sortTags < theTags.getTagCount(); sortTags++) {//fill orderedMasterTags with indices (index 1) and read counts (index 0)
           orderByRead[sortTags][0]= theTags.getReadCount(sortTags);
           orderByRead[sortTags][1]= sortTags;
       }

       Arrays.sort(orderByRead, new java.util.Comparator<int[]>() {//sort 2d array by read count, descending
           @Override
           public int compare(int[] o1, int[] o2) {
               if(o1[0]>o2[0]) return -1; // -1 for descending order
               else if(o1[0]<o2[0]) return 1; // 1 for descending order
               else return 0;
           }
       });
       
       int[][] stats= new int[(theTags.getTagCount()/divisor)+1][23];//index 0 (which division), 1 (alignCount), 2(multiMap), 3(telomeric), 4(pericentromeric), 5(centromeric, 6-18(chromosome count)
       int[][] chrPos= {//where chr1 is index 0, chr2 is index 1 and so on 
               {133300000,89300000,94600000,104200000,101600000,49800000,55300000,45900000,68600000,59300000},//holds the starting centromeric position for each chromosome using the map position of functional centromeres from Wolfgruber at al Table 1 
               {133900000,91100000,95400000,105000000,108600000,50400000,55700000,48000000,69200000,60700000},//holds the ending centromeric position for each chromosome using the map position of functional centromeres from Wolfgruber at al Table 1. Use the greatest range for chr5
               {100000000,65000000,30000000,40000000,40000000,35000000,50000000,25000000,35000000,25000000},//holds starting "pericentromeric" position extrapolated from Gore et al 2009 based on where the recombination rate decreases
               {150000000,160000000,120000000,130000000,140000000,85000000,90000000,100000000,85000000,75000000}//holds ending "pericentromeric" position extrapolated from Gore et al 2009 based on where the recombination rate increases
           };
       
       int division= 0;
       
       int realTBTIndex;
       long[] tag;
       int TOPMIndex;
       int chr;
       int minPos;
       
       for (int orderedTBTIndex= 0; orderedTBTIndex < theTags.getTagCount(); orderedTBTIndex++) {
           realTBTIndex= orderByRead[orderedTBTIndex][1];
           tag= theTags.getTag(realTBTIndex);
           TOPMIndex= theTOPM.getTagIndex(tag);
           chr=  theTOPM.getChromosome(TOPMIndex);
           minPos= theTOPM.getStartPosition(TOPMIndex);
                      
           if (theTOPM.getStartPosition(TOPMIndex) > -1) {
               stats[division][chr+6]++;//plus 6 because 6 other variables held in array (for chromosomes 0-12)
               stats[division][1]++;//num of tags that align
               
               if (chr > 0 && chr < 11) {
                   if ((minPos < chrPos[1][chr-1]) && (minPos > chrPos[0][chr-1])) stats[division][5]++;                     
                   else if (minPos < chrPos[3][chr-1] && minPos > chrPos[2][chr-1]) stats[division][4]++;                     
                   else stats[division][3]++;
               }
           }
           
           if (theTOPM.getMultiMaps(TOPMIndex) > 1) stats[division][2]++;
           if (theTOPM.getMultiMaps(TOPMIndex) == 1) stats[division][19]++;
           if (theTOPM.getMultiMaps(TOPMIndex) == 2) stats[division][20]++;
           if (theTOPM.getMultiMaps(TOPMIndex) > 2 && theTOPM.getMultiMaps(TOPMIndex) != 99) stats[division][21]++;
           if (theTOPM.getMultiMaps(TOPMIndex) == 99) stats[division][22]++;
           
           if (((orderedTBTIndex)%divisor == 0 && orderedTBTIndex != 0) || theTags.getTagCount()-1 == orderedTBTIndex) {
               division++;
               stats[division-1][0]= division;
               System.out.println("finish "+(division*divisor)+" tags in the TBT, ordered by read count");
               
           }
       }
       try {
                    outputFile.createNewFile();
                    outputFile.setWritable(true);
                    BufferedWriter BW = new BufferedWriter(new FileWriter(outputFile));

                    BW.write("Division\tAligned\tMultipleMapping\tNumAlignedTelomeric\tNumAlignedPericentromeric\tNumAlignedCentromeric");
                    for (int i= 0; i < 13; i++) {
                        BW.write("\tChr"+i);
                    }
                    BW.write("\tSingleMapping\tTwoPositions\tGreaterThanTwoPositions\tUnknownNumMultipleMappings\n");
                    
                    for (int div=0; div < stats.length; div++) {
                        for (int content= 0; content < 23; content++) {
                            BW.write(stats[div][content]+"\t");                            
                        }
                        BW.write("\n");
                    }
                    BW.close();
                }

                catch (IOException e) {
                    System.out.println(e);
                }
   }
      
       public static void HapmapToCHIAMO(String inFile) { //the genotypic file for use in IMPUTE2 and Beagle. Missing data is given the HW probability based on MajorMinor Allele Freq
        String hapMapFileName= dir+inFile+".hmp.txt";
        String outCHIAMOFile= dir+inFile+".gens";
        Alignment a= ImportUtils.readFromHapmap(hapMapFileName, null);
        DecimalFormat df= new DecimalFormat("0.####");
        
        try{
          DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outCHIAMOFile), 655360));
          
          for (int site= 0; site < a.getSiteCount(); site++) {
              outStream.writeBytes(a.getLocusName(site) +" "+a.getSNPID(site)+" "+a.getPositionInLocus(site)+" "+a.getMajorAlleleAsString(site)+" "+a.getMinorAlleleAsString(site));
              double p= a.getMajorAlleleFrequency(site);
              double q= a.getMinorAlleleFrequency(site);
              double f;
              
              for (int taxa= 0; taxa<a.getSequenceCount(); taxa++) {
                  f= ((double) a.getHeterozygousCountForTaxon(taxa)/((double) a.getSiteCount()-((double) a.getTotalGametesNotMissingForTaxon(taxa)/2))>.015)?.9:.2;
                  if (a.isHeterozygous(taxa, site) == true) outStream.writeBytes(" 0 1 0");
                  else if (a.getBaseArray(taxa, site)[0] == a.getMajorAllele(site) && a.getBaseArray(taxa, site)[1] == a.getMajorAllele(site)) outStream.writeBytes(" 1 0 0");
                  else if (a.getBaseArray(taxa, site)[0] == a.getMajorAllele(site) || a.getBaseArray(taxa, site)[1] == a.getMajorAllele(site)) outStream.writeBytes(" "+df.format(p*f)+" "+df.format(q*(1-f))+" 0");
                  else if (a.getBaseArray(taxa, site)[0] == a.getMinorAllele(site) && a.getBaseArray(taxa, site)[1] == a.getMinorAllele(site)) outStream.writeBytes(" 0 0 1");
                  else if (a.getBaseArray(taxa, site)[0] == a.getMinorAllele(site) || a.getBaseArray(taxa, site)[1] == a.getMinorAllele(site)) outStream.writeBytes(" 0"+df.format(p*(1-f))+" "+df.format(q*f));
                  else outStream.writeBytes(" "+df.format(Math.pow(p,2)*f)+" "+df.format(2*p*q*(1-f))+" "+df.format(Math.pow(q, 2)*f));
              }
              outStream.writeBytes("\n");
           }
          outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
        
    }
       
   public static void HapmapToSample(String inFile) { //the basic sample file for use in IMPUTE2. Calculates missing based on diploid value
        String hapMapFileName= dir+inFile+".hmp.txt";
        String outSamplesFile= dir+inFile+".samples";
        Alignment a= ImportUtils.readFromHapmap(hapMapFileName, null);
        DecimalFormat df= new DecimalFormat("0.##");
        double miss;
        
        try{
          DataOutputStream outStream= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(outSamplesFile), 655360));
          outStream.writeBytes("ID_1 ID_2 missing\n0 0 0");
          
          for (int taxon= 0; taxon < a.getSequenceCount(); taxon++) {
              miss= 1.0-((double) a.getTotalGametesNotMissingForTaxon(taxon)/((double) a.getSiteCount()*2.0));
              outStream.writeBytes("\n"+a.getTaxaName(taxon)+" "+a.getTaxaName(taxon)+" "+df.format(miss));              
           }
          outStream.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
     }
   
   public static void LegendAndHaps(Alignment a) { //to generate a halotype and accompanying legends file from an alignment for IMPUTE2. Throws out any alignments with hets
       String hapsFileName= dir+a+".haps";
       String legendFileName= dir+a+".legend";
       
       try{
          DataOutputStream outStreamHaps= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(hapsFileName), 655360));
          DataOutputStream outStreamLegend= new DataOutputStream(new BufferedOutputStream(new FileOutputStream(legendFileName), 655360));
          outStreamLegend.writeBytes("rdID position a0 a1\n");
          
          for (int site= 0; site < a.getSiteCount(); site++) {
              
           }
          outStreamHaps.close();
          outStreamLegend.close();
       }
        
      catch(IOException e) {
           System.out.println(e);
       }
   }
   
   public static void SubsetHapmapByPosition(String inFile, int startPos, int endPos, boolean gz) {
        String inHapMapFileName= (gz==true)?dir+inFile+".hmp.txt.serial.gz":dir+inFile+".hmp.txt";
        String outHapMapFileName= dir+inFile+"subset_"+startPos+"-"+endPos+".hmp.txt";
        Alignment a= (gz==true)?ImportUtils.readAlignmentFromSerialGZ(inHapMapFileName):ImportUtils.readFromHapmap(inHapMapFileName, null);
        MutableNucleotideAlignment newAlign= MutableNucleotideAlignment.getInstance(a);
        for (int site= 0; site<a.getSiteCount(); site++) {
            if (a.getPositionInLocus(site)<startPos || a.getPositionInLocus(site)>endPos) newAlign.clearSiteForRemoval(site);
        }
        newAlign.clean();
        System.out.println("First site pos: "+newAlign.getPositionInLocus(0)+"\nLast site pos: "+newAlign.getPositionInLocus(newAlign.getSiteCount()-1));
        ExportUtils.writeToHapmap(newAlign, true, outHapMapFileName, '\t', null);
   }
   
   public static Alignment SubsetHapmapByMAF(Alignment a, double MAF, boolean toFile) {
       System.out.println("Subsetting by MAF "+MAF);
       System.out.println("Original number of sites: "+a.getSiteCount());
       MutableNucleotideAlignment newAlign= MutableNucleotideAlignment.getInstance(a);
        for (int site= 0; site<a.getSiteCount(); site++) {
            if (a.getMinorAlleleFrequency(site)<MAF) newAlign.clearSiteForRemoval(site);
        }
        newAlign.clean();
        System.out.println("Sites remaining: "+newAlign.getSiteCount());
        if (toFile==true) ExportUtils.writeToHapmap(newAlign, true, dir+a+"subsetByMAF_"+MAF+".hmp.txt.gz", '\t', null);
        return newAlign;
   }
   
   public static Alignment FilterForPolymorphicSites(Alignment a) {
       System.out.println("Filter for polymorphic sites...");
       ArrayList<Integer> subSite= new ArrayList<Integer>();
        for (int s= 0; s<a.getSiteCount(); s++) {
            if (a.isPolymorphic(s) == true) subSite.add(s);
        }
        int[] keepSite= ArrayUtils.toPrimitive(subSite.toArray(new Integer[NumPolymorphicSites(a)]));
        Alignment newAlign= FilterAlignment.getInstance(a, keepSite);
        System.out.println("Final num of taxa: "+newAlign.getSequenceCount());
        System.out.println("Final num of sites: "+newAlign.getSiteCount());
        return newAlign;
   }
   
   public static int NumPolymorphicSites(Alignment a) {
        int poly= 0;
        for (int s= 0; s<a.getSiteCount(); s++) {
            if (a.isPolymorphic(s) == true) poly++;
        }
        return poly;
   }
   
   public static Alignment SubsetHapmapByTaxaCov(Alignment a, double minCov, boolean toFile, boolean filterPolymorphic) {
        System.out.println("Subsetting by minimum taxa coverage: "+minCov);
        System.out.println("Original num of taxa: "+a.getSequenceCount());
        System.out.println("Original num of sites: "+a.getSiteCount());
        String outHapMapFileName= dir+a+"subset_"+"_minCov"+minCov+".hmp.txt.gz";
        IdGroup IDs= a.getIdGroup();
        boolean[] badTaxa= new boolean[a.getSequenceCount()];
        double numSites= a.getSiteCount();
        for (int taxon= 0; taxon<a.getSequenceCount(); taxon++) {
            if ((((double) a.getTotalNotMissingForTaxon(taxon))/numSites) <= minCov) badTaxa[taxon]= true;
        }
        IdGroup removeIDs= IdGroupUtils.idGroupSubset(IDs, badTaxa);
        Alignment align= FilterAlignment.getInstanceRemoveIDs(a, removeIDs);
        if (filterPolymorphic==true) {
            //filter for sites that are polymorphic
            Alignment newAlign= FilterForPolymorphicSites(align);
            if (toFile==true) ExportUtils.writeToHapmap(newAlign, true, outHapMapFileName, '\t', null);
            return newAlign;
        }
        else {
            if (toFile==true) ExportUtils.writeToHapmap(align, true, outHapMapFileName, '\t', null);
            System.out.println("Final num of taxa: "+align.getSequenceCount());
            System.out.println("Final num of sites: "+align.getSiteCount());
            return align;
        }
   }
    public static String[] readInTxtNames(String inFile, boolean permissive) {
       Set<String> names= new HashSet<String>();
       try {
            FileInputStream fis= new FileInputStream(inFile);
            Scanner scanner= new Scanner(fis);
            do {
                String next= scanner.nextLine();
                if (permissive==true) names.add(next.substring(0, next.indexOf(":")));
                else names.add(next);
            }
            while (scanner.hasNextLine());
            scanner.close();
            fis.close();
        }
        catch (Exception e) {
            System.out.println("Problem reading in taxa names");
        }
        ArrayList<String> sortNames= new ArrayList<String>();
        sortNames.addAll(names);
        Collections.sort(sortNames);
        String[] nameArray= Arrays.copyOf(sortNames.toArray(),sortNames.size(),String[].class);
        return nameArray;
   }
   
   public static void MergeAlignments(String inFile1, boolean gz1, String inFile2, boolean gz2, String outFileName) {
       Alignment a1= ImportUtils.readFromHapmap(dir+inFile1+(gz1==true?".hmp.txt.gz":".hmp.txt"), null);
       Alignment a2= ImportUtils.readFromHapmap(dir+inFile2+(gz2==true?".hmp.txt.gz":".hmp.txt"), null);
       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a1, a1.getSequenceCount()+a2.getSequenceCount(), a1.getSiteCount());
       IdGroup a2ID= a2.getIdGroup();
       for (int taxon = 0; taxon < a2.getSequenceCount(); taxon++) {
           mna.addTaxon(a2ID.getIdentifier(taxon));
           for (int site = 0; site < a1.getSiteCount(); site++) {
               mna.setBase(a1.getSequenceCount()+taxon, site, a2.getBase(taxon, site));
           }
       }
       mna.clean();
       ExportUtils.writeToHapmap(mna, true, outFileName+".hmp.txt.gz", '\t', null);
   }
   
   public static Alignment SubsetHapmapByHeterozygosity(String inHapmap, boolean gz, double hetCutoff,boolean desireHets, boolean toFile, boolean filterPolymorphic) {
        Alignment a= ImportUtils.readFromHapmap(dir+inHapmap+(gz==true?".hmp.txt.gz":".hmp.txt"), null);
        System.out.println("Subsetting by Heterozygosity "+(desireHets==true?"greater than "+hetCutoff:"less than or equal to "+hetCutoff));
        System.out.println("Original num of taxa: "+a.getSequenceCount());
        System.out.println("Original num of sites: "+a.getSiteCount());
        String outHapMapFileName= (desireHets==false)?dir+inHapmap+"subset_lessThanEqual"+hetCutoff+"Het.hmp.txt.gz":dir+inHapmap+"subset_greaterThan"+hetCutoff+"Het.hmp.txt.gz";
        
        IdGroup IDs= a.getIdGroup();
        boolean[] badTaxa= new boolean[a.getSequenceCount()];
        double numSites= a.getSiteCount();
        if (desireHets==false) {
            for (int taxon= 0; taxon<a.getSequenceCount(); taxon++) {
            if ((((double) a.getHeterozygousCountForTaxon(taxon))/numSites) > hetCutoff) badTaxa[taxon]= true;
            }
        }
        else {
            for (int taxon= 0; taxon<a.getSequenceCount(); taxon++) {
            if ((((double) a.getHeterozygousCountForTaxon(taxon))/numSites) <= hetCutoff) badTaxa[taxon]= true;
            }
        }
        IdGroup removeIDs= IdGroupUtils.idGroupSubset(IDs, badTaxa);
        Alignment align= FilterAlignment.getInstanceRemoveIDs(a, removeIDs);
        if (filterPolymorphic==true) {
            //filter for sites that are polymorphic 
            Alignment newAlign= FilterForPolymorphicSites(align);
            if (toFile==true) ExportUtils.writeToHapmap(newAlign, true, outHapMapFileName, '\t', null);
            return newAlign;
        }
        else {
            if (toFile==true) ExportUtils.writeToHapmap(align, true, outHapMapFileName, '\t', null);
            System.out.println("Final num of taxa: "+align.getSequenceCount());
            System.out.println("Final num of sites: "+align.getSiteCount());
            return align;
        }
   }   
   
   public static void RemoveSitesWithGaps(String inFile, boolean gz) {
       String inHapMapFileName= (gz==true)?dir+inFile+".hmp.txt.gz":dir+inFile+".hmp.txt";
       Alignment a= ImportUtils.readFromHapmap(inHapMapFileName, null);
       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
       for (int site = 0; site < a.getSiteCount(); site++) {
           if (mna.getMajorAllele(site)==NucleotideAlignmentConstants.GAP_ALLELE||mna.getMajorAllele(site)==NucleotideAlignmentConstants.INSERT_ALLELE
               ||mna.getMinorAllele(site)==NucleotideAlignmentConstants.GAP_ALLELE||mna.getMinorAllele(site)==NucleotideAlignmentConstants.INSERT_ALLELE)
               mna.clearSiteForRemoval(site);
       }
       mna.clean();
       ExportUtils.writeToHapmap(mna, true, dir+inFile+"IndelsRemoved.hmp.txt.gz", '\t', null);
   }
   
   //FIX THIS TO CONFORM TO NEW LD CODE 7/16/2013
//   public static void CheckSitesForLD(String posToCheckFile, boolean posToCheckGZ, String inFile, boolean inGz, double r2cutoff) {
//       String posFileName= (posToCheckGZ==true)?dir+posToCheckFile+".hmp.txt.gz":dir+posToCheckFile+".hmp.txt";
//       Alignment siteAlign= ImportUtils.readFromHapmap(posFileName, null);
//       String inHapMapFileName= (inGz==true)?dir+inFile+".hmp.txt.gz":dir+inFile+".hmp.txt";
//       Alignment a= ImportUtils.readFromHapmap(inHapMapFileName, null);
//       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(siteAlign);
//       int[] pos= siteAlign.getPhysicalPositions();
//       for (int sites:pos) {
//           int focus= a.getSiteOfPhysicalPosition(sites, null);
//           LinkageDisequilibrium ld= new LinkageDisequilibrium(a, 10,
//                   LinkageDisequilibrium.testDesign.SlidingWindow, focus, null, false, -1, null);
//           ld.run();
//           double[] r2= new double[18];
//           for (int site = 1; site < 10; site++) {
//               r2[site]= ld.getRSqr(focus, focus-site);
//               r2[r2.length-site]= ld.getRSqr(focus, focus+site);
//           }
//           for(double r:r2){System.out.println(r);}
//           Arrays.sort(r2);
//           if (r2[10]>r2cutoff) mna.clearSiteForRemoval(sites);
//       }
//       mna.clean();
//       ExportUtils.writeToHapmap(mna, true, posToCheckFile+"_goodLD.hmp.txt.gz", '\t', null);
//   }
   
   //takes taxa from inFile and sites from refFile. Designed to compare 55k
   public static void CheckSitesForIdentityByPosition(String inFile, boolean gz, String refFile, boolean refGz, double siteThreshold, double taxonThreshold) {
       String inHapMapFileName= (gz==true)?dir+inFile+".hmp.txt.gz":dir+inFile+".hmp.txt";
       Alignment a= ImportUtils.readFromHapmap(inHapMapFileName, null);
       
       String inRefFileName= (refGz==true)?dir+refFile+".hmp.txt.gz":dir+refFile+".hmp.txt";
       Alignment ref= ImportUtils.readFromHapmap(inRefFileName, null);
       MutableNucleotideAlignment mnaRef= MutableNucleotideAlignment.getInstance(ref);
       
       int[] refPos= ref.getPhysicalPositions();
       int[] matchTaxon= new int[a.getSequenceCount()];//holds the index of corresponding taxon in ref
       String[] refNames= new String[ref.getSequenceCount()];
       for (int taxon = 0; taxon < refNames.length; taxon++) {
           refNames[taxon]= ref.getIdGroup().getIdentifier(taxon).getNameLevel(0);
       }      
       for (int taxon = 0; taxon < matchTaxon.length; taxon++) {
           String unkName= a.getIdGroup().getIdentifier(taxon).getNameLevel(0);
           matchTaxon[taxon]= Arrays.binarySearch(refNames, unkName);
       }
       double[] badSites= new double[refPos.length];
       double[] halfSites= new double[refPos.length];
       double[] testedSites= new double[refPos.length];
       double[] badTaxa= new double[a.getSequenceCount()];
       double[] halfTaxa= new double[a.getSequenceCount()];
       double[] testedTaxa= new double[a.getSequenceCount()];
       for (int site= 0;site<a.getSiteCount();site++) {
           int refIndex= Arrays.binarySearch(refPos, a.getPositionInLocus(site));
           if (refIndex<0) continue;
           int refSite= ref.getSiteOfPhysicalPosition(a.getPositionInLocus(site), null);
           for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
               byte inBase= a.getBase(taxon, site);
               byte refBase= ref.getBase(matchTaxon[taxon], refSite);
               byte[] inBaseArray= a.getBaseArray(taxon, site);
               byte[] refBaseArray= ref.getBaseArray(matchTaxon[taxon], refSite);
               if (inBase==diploidN||refBase==diploidN) continue;
               testedSites[refIndex]+=1.0;
               testedTaxa[taxon]+=1.0;
               if (AlignmentUtils.isEqual(inBaseArray, refBaseArray)==true) continue;
               if (a.isHeterozygous(taxon, site)==true||ref.isHeterozygous(matchTaxon[taxon], refSite)==true) {
                   if (inBaseArray[0]==refBaseArray[0]||inBaseArray[1]==refBaseArray[0]||inBaseArray[0]==refBaseArray[1]||
                        inBaseArray[1]==refBaseArray[1]) {
                       halfSites[refIndex]+=1.0;
                       halfTaxa[taxon]+=1.0;
                   }
               }
               else {
                   badSites[refIndex]+=1.0;
                   badTaxa[taxon]+=1.0;
               }
           }
           if (badSites[refIndex]/testedSites[refIndex]>siteThreshold) mnaRef.clearSiteForRemoval(refIndex);
           System.out.println(a.getSNPID(site)+"\t"+badSites[refIndex]/testedSites[refIndex]
                   +"\t"+halfSites[refIndex]/testedSites[refIndex]);
       }
       boolean[] badIDs= new boolean[a.getSequenceCount()];
       for (int taxon = 0; taxon < a.getSequenceCount(); taxon++) {
           System.out.println(a.getTaxaName(taxon)+"\t"+badTaxa[taxon]/testedTaxa[taxon]
                   +"\t"+halfTaxa[taxon]/testedTaxa[taxon]);
           if (badTaxa[taxon]/testedTaxa[taxon]>taxonThreshold) badIDs[taxon]= true;
           
       }
       IdGroup removeIDs= IdGroupUtils.idGroupSubset(a.getIdGroup(), badIDs);
       Alignment goodTaxa= FilterAlignment.getInstanceRemoveIDs(a, removeIDs);
       System.out.println("Write "+goodTaxa.getSequenceCount()+" taxa to file");
       ExportUtils.writeToHapmap(goodTaxa, true, dir+inFile+"_cleanTaxa.hmp.txt.gz", '\t', null);
       
       mnaRef.clean();
       System.out.println("Write "+mnaRef.getSiteCount()+" sites to file");
       ExportUtils.writeToHapmap(mnaRef, true, dir+refFile+"_cleanSites.hmp.txt.gz", '\t', null);
   }
   public static void SitesWithSameNamesOrPositions(String inFile, boolean gz) {
       String inHapMapFileName= (gz==true)?dir+inFile+".hmp.txt.gz":dir+inFile+".hmp.txt";
       Alignment a= ImportUtils.readFromHapmap(inHapMapFileName, null);
       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
       String[] SNPs= a.getSNPIDs();
       int[] pos= a.getPhysicalPositions();
       System.out.println("duplicate sites set for removal:");
       for (int id1 = 0; id1 < SNPs.length; id1++) {
           for (int id2 = id1+1; id2 < SNPs.length; id2++) {
               if (SNPs[id1]==SNPs[id2]||pos[id1]==pos[id2]) {
                   mna.clearSiteForRemoval(id2);
                   System.out.println(SNPs[id2]+" ("+pos[id2]+")"+" ("+id2+" matches "+id1+")");
               }
           }
       }
       mna.clean();
       ExportUtils.writeToHapmap(mna, true, dir+inFile+"DupSitesRemoved.hmp.txt.gz", '\t', null);
   }
   
   public static void MaskSites(String inFile, boolean gz, int maskRate) {
       String inHapMapFileName= (gz==true)?dir+inFile+".hmp.txt.gz":dir+inFile+".hmp.txt";
       String outHapMapFileName= dir+inFile+"_masked.hmp.txt";
       Alignment a= ImportUtils.readFromHapmap(inHapMapFileName, null);
       System.out.println("Read in hapmap "+inHapMapFileName);
       MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
       for (int taxon= 0;taxon<a.getSequenceCount();taxon++) {            
            for (int site= taxon;site<a.getSiteCount();site+= maskRate) {
                mna.setBase(taxon, site, diploidN);
            }
       }
       mna.clean();
       ExportUtils.writeToHapmap(mna, true, outHapMapFileName, '\t', null);
   }
   
   public static void HapmapToVCF(String inFile, boolean gz, boolean removeGapSites) {
       Alignment a= ImportUtils.readFromHapmap(dir+inFile+(gz==true?".hmp.txt.gz":".hmp.txt"),null);
       if (removeGapSites==true) {
            System.out.println("Sites in original alignment: "+a.getSiteCount());
            MutableNucleotideAlignment mna= MutableNucleotideAlignment.getInstance(a);
            for (int site= 0; site<a.getSiteCount();site++) {
                if (a.getAlleles(site)[2]==NucleotideAlignmentConstants.GAP_DIPLOID_ALLELE) mna.clearSiteForRemoval(site);
            }
            mna.clean();
            System.out.println("Sites after filtering: "+mna.getSiteCount());
       }
       ExportUtils.writeToVCF(a, dir+inFile+".vcf", '\t');
   }
   
   public static void RandomlySampleAlignment(String inFile, boolean gz, int numTaxa, boolean filterPolymorphic) {
       String outHapMapFileName= dir+inFile+"RndSample"+numTaxa+".hmp.txt.gz";
       Alignment a= ImportUtils.readFromHapmap(dir+inFile+(gz==true?".hmp.txt.gz":".hmp.txt"),null);
       IdGroup IDs= a.getIdGroup();
        boolean[] badTaxa= new boolean[a.getSequenceCount()];
        int count= 0;
        while (count<=a.getSequenceCount()-numTaxa-1) {
            int num= (int)Math.floor(Math.random()*a.getSequenceCount());
            if (badTaxa[num]==false) count++;
            badTaxa[num]= true;
        }
        IdGroup removeIDs= IdGroupUtils.idGroupSubset(IDs, badTaxa);
        Alignment align= FilterAlignment.getInstanceRemoveIDs(a, removeIDs);
        if (filterPolymorphic==true) {
            //filter for sites that are polymorphic 
            align= FilterForPolymorphicSites(align);
        }
        ExportUtils.writeToHapmap(align, true, outHapMapFileName, '\t', null);
   }
   
   public static void SubsetByTaxaName(String forSubset, String forIDGroup) {
       Alignment toSubset= ImportUtils.readFromHapmap(dir+forSubset+".hmp.txt.gz", null);
       Alignment useIDGroup= ImportUtils.readFromHapmap(dir+forIDGroup+".hmp.txt.gz", null);
       String outHapMapFileName= dir+forSubset+"subsetBy"+forIDGroup+".hmp.txt.gz";
       IdGroup match= useIDGroup.getIdGroup();
       IdGroup sub= toSubset.getIdGroup();
       boolean[] badTaxa= new boolean[toSubset.getSequenceCount()];
       String[] IDString= new String[useIDGroup.getSequenceCount()];
       
       for (int taxon = 0; taxon < useIDGroup.getSequenceCount(); taxon++) {
           IDString[taxon]= match.getIdentifier(taxon).getNameLevel(0);
       }
       for (int taxon = 0; taxon < toSubset.getSequenceCount(); taxon++) {
           String name= sub.getIdentifier(taxon).getNameLevel(0);
           int index= Arrays.binarySearch(IDString, name);
           if (index<0) badTaxa[taxon]= true;
       }
       IdGroup removeIDs= IdGroupUtils.idGroupSubset(sub, badTaxa);
       Alignment align= FilterAlignment.getInstanceRemoveIDs(toSubset, removeIDs);
       ExportUtils.writeToHapmap(align, true, outHapMapFileName, '\t', null);
   }
   
   public static void subsetHDF5FromTxt(String inFileRoot, String taxaNamesRoot, boolean permissive, boolean h5) {
       MutableNucleotideAlignmentHDF5 mnah5= (MutableNucleotideAlignmentHDF5)ImportUtils.readGuessFormat(dir+inFileRoot+".hmp.h5", false);
       Set<String> names= new HashSet<String>();
       IdGroup h5IDs= mnah5.getIdGroup();
       if (permissive==true) System.out.println("permissive mode on (matches only the first part of the name before the first colon - will select duplicate samples of taxa)");
        try {
            FileInputStream fis= new FileInputStream(dir+taxaNamesRoot+".txt");
            Scanner scanner= new Scanner(fis);
            while (scanner.hasNextLine()) {
                String next= scanner.nextLine();
                if (permissive==true) names.add(next.substring(0, next.indexOf(":")));
                else names.add(next);
            }
            scanner.close();
            fis.close();
        }
        catch (Exception e) {
        }
        ArrayList<String> sortNames= new ArrayList<String>();
        sortNames.addAll(names);
        Collections.sort(sortNames);
        String[] nameArray= Arrays.copyOf(sortNames.toArray(),sortNames.size(),String[].class);
        ArrayList<Identifier> keep= new ArrayList<Identifier>();
        for (int h5Taxa = 0; h5Taxa < h5IDs.getIdCount(); h5Taxa++) {
            String currTaxon= permissive==true?h5IDs.getIdentifier(h5Taxa).getName():h5IDs.getIdentifier(h5Taxa).getFullName();
            if (Arrays.binarySearch(nameArray, currTaxon)>-1) keep.add(h5IDs.getIdentifier(h5Taxa));
        }
        SimpleIdGroup subIDs= new SimpleIdGroup(keep);
        Alignment subset= FilterAlignment.getInstance(mnah5, subIDs);
        System.out.println("subsetting "+subset.getSequenceCount()+" taxa from "+h5IDs.getIdCount()+" taxa based on "+names.size()+" unique names");
        if (h5==true) {
            if (permissive==true) ExportUtils.writeToHDF5(subset, dir+inFileRoot+"PermissiveSubsetBy"+taxaNamesRoot+".hmp.h5");
            else ExportUtils.writeToHDF5(subset, dir+inFileRoot+"StrictSubsetBy"+taxaNamesRoot+".hmp.h5");
            System.out.println("writing to hdf5");
        }
        else {
            if (permissive==true) ExportUtils.writeToHapmap(subset, true, dir+inFileRoot+"PermissiveSubsetBy"+taxaNamesRoot+".hmp.h5", '\t', null);
            else ExportUtils.writeToHapmap(subset, true, dir+inFileRoot+"StrictSubsetBy"+taxaNamesRoot+".hmp.h5", '\t', null);
            System.out.println("writing to hapmap");
        }
   }
   
   //from Alberto to make a beagle 3 file
      public static String saveDelimitedAlignment(Alignment theAlignment, String delimit, String saveFile) {

        if ((saveFile == null) || (saveFile.length() == 0)) {
            return null;
        }
        saveFile = Utils.addSuffixIfNeeded(saveFile, ".txt");
        FileWriter fw = null;
        BufferedWriter bw = null;
        try {

            fw = new FileWriter(new File(saveFile));
            bw = new BufferedWriter(fw);

            bw.write("I\tid");
            int numSites = theAlignment.getSiteCount();
            for (int j = 0; j < theAlignment.getSequenceCount(); j++) {
                bw.write(delimit);
                bw.write(theAlignment.getIdGroup().getIdentifier(j).getFullName());
                bw.write(delimit);
                bw.write(theAlignment.getIdGroup().getIdentifier(j).getFullName());              
            }
            bw.write("\n");
            for (int i = 0; i < numSites; i++) {
                bw.write("M\t");
//                bw.write(String.valueOf(theAlignment.getPositionInLocus(i)));
                    bw.write(theAlignment.getSNPID(i));
                for (int r = 0, n = theAlignment.getSequenceCount(); r < n; r++) {
                    bw.write(delimit);
                    bw.write(theAlignment.getBaseAsStringArray(r, i)[0]);
                    bw.write(delimit);
                    bw.write(theAlignment.getBaseAsStringArray(r, i)[1]);
                }
                bw.write("\n");
            }

            return saveFile;

        } catch (Exception e) {
            throw new IllegalArgumentException("Error writing Delimited Alignment: " + saveFile + ": " + ExceptionUtils.getExceptionCauses(e));
        } finally {
            try {
                bw.close();
                fw.close();
            } catch (Exception e) {
                // do nothing
            }
        }

    }
      
    //combine alignments, meant for the seqToGenos_partX.hmp.h5 syntax
   public static void combineAlignmentsKeepDepth(String inFileRoot, int startNum, int endNum) {
       MutableNucleotideAlignmentHDF5[] aligns= new MutableNucleotideAlignmentHDF5[endNum+1-startNum];
       for (int i = startNum; i < endNum+1; i++) {aligns[i]= MutableNucleotideAlignmentHDF5.getInstance(inFileRoot+i+"hmp.h5");}
       Alignment newAlign= MutableNucleotideDepthAlignment.getInstance(aligns);
       ExportUtils.writeToMutableHDF5(newAlign, inFileRoot+"_combined"+startNum+"-"+endNum+".hmp.h5", null, true);
   }

   
   public static void main (String args[]) {
       TasselPrefs.putAlignmentRetainRareAlleles(false);
             
//       //args for FilterMergedTBT
//       int minSharedTaxa= 2;
//       String inTBTFileName= dir+ "/landraceCombo/mergedTBT/SW_RI_Span_landrace_min1.tbt.byte";//input tbt file pathway
//       String outTBTFileName= dir+ "/landraceCombo/mergedTBT/SW_RI_Span_landrace_min"+minSharedTaxa+".tbt.byte";//output tbt file
//       FilterMergedTBT(minSharedTaxa, inTBTFileName, outTBTFileName);
       
       //args for CharTOPM
//       String inTOPMFileName= dir+"/landraceCombo/topm/SW_RI_Span_landraces_min2_bowtie2MissingOrganelles.topm";
//       String inTBTFileName= dir+"/landraceCombo/mergedTBT/SW_RI_Span_landrace_min1.tbt.byte";
//       int divisor= 1000000;
//       CharTOPMWithTBT(inTOPMFileName, inTBTFileName, divisor);   
       
//       //args for CharTOPMWithTagCount
//       String inTOPMFileName= dir+"/AllZeaMasterTags_c10_20120703.topm";
//       String inTagCountFileName= dir+"/AllZeaMasterTags_c10_20120606.cnt";
//       int divisor= 1000000;
//       CharTOPMWithTagCount(inTOPMFileName, inTagCountFileName, divisor);
       
       //args for subsetting and formatting for Impute2
//       int startPos= 30000000;
//       int endPos= 40000000;
//       dir= "/Users/kelly/Documents/GBS/FinalRev1_BPECFilteredSNPsSubset/";
//       String inGBSFile= "AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_Rev1_chr8";
//       SubsetHapmapByPosition(inGBSFile,startPos,endPos,true);
//       inGBSFile= "AllTaxa_BPEC_AllZea_GBS_Build_July_2012_FINAL_Rev1_chr8subset_"+startPos+"-"+endPos;
//       HapmapToCHIAMO(inGBSFile);
//       HapmapToSample(inGBSFile);
//       
//       dir= "/home/local/MAIZE/kls283/GBS/Imputation/minCov10percentsubsets04_PivotMergedTaxaTBT.c10/";
//       String inHapmap= "04_PivotMergedTaxaTBT.c10_s0_s4095subset__minCov0.1_HaplotypeMerge";
//       SubsetHapmapByTaxaCov(inHapmap,0.6,false);
//       SubsetHapmapByPosition(inWGSFile,startPos,endPos,false);
//       inWGSFile= "maizeHapMapV2_B73RefGenV2_201203028_chr8subset"+startPos+"-"+endPos;
//       HapmapToCHIAMO(inWGSFile);
//       HapmapToSample(inWGSFile);
       
//       dir= "/home/local/MAIZE/kls283/GBS/Imputation/2.6_MergeDupSNPs/";
//       String inHapmap= "AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1";
////       SubsetHapmapByTaxaCov(inHapmap, .1, true);
//       MaskSites(inHapmap,false,300);
       
       //subset for heterozygosity
       dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
//       String inHapmap= "AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1";
//       SubsetHapmapByHeterozygosity(inHapmap,true,.01,true,true,false);
       
//       dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
//       String inFile= "04_PivotMergedTaxaTBT.c10_s0_s24575subset__minCov0.1";
//       SitesWithSamePhysicalPositions(inFile,false);
       
//       dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
//       String inFile= "AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_minHet0.024";
//       HapmapToVCF(inFile, true, true);
       
//       Alignment a= ImportUtils.readFromHapmap(dir+"maizeHapMapV2_B73RefGenV2_201203028_chr10.hmp.txt", null);
//       saveDelimitedAlignment(a,"\t",dir+"maizeHapMapV2_B73RefGenV2_201203028_chr10.bgl");
       
//       RandomlySampleAlignment("SEED_12S_GBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1RndSample1000",true,500,false);
       
//       String toSubset= "AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1";
//       String IDGroup= "SNP55K_maize282_AGPv2_20100513_1.chr10";
//       SubsetByTaxaName(toSubset, IDGroup);
       
//       dir= "/Users/kelly/Documents/GBS/Imputation/SmallFiles/";
//       SitesWithSameNamesOrPositions("SNP55K_maize282_AGPv2_20100513_1",false);
//       
//       dir= "/Users/kelly/Documents/GBS/Imputation/SmallFiles/";
//       CheckSitesForLD("RIMMA_282_SNP55K_AGPv2_20100513__S45391.chr10_matchTo_RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1",true,
//               "RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.2",true, .5);
       
//       CheckSitesForIdentityByPosition("RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.2",true,"RIMMA_282_SNP55K_AGPv2_20100513__S45391.chr10_matchTo_RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1",true,.1,.1);
//       dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
//       subsetHDF5FromTxt("AllZeaGBSv27i3b.imp","Ames(no EP or GEM)",false,false);
//       
//       dir= "/Users/kelly/Documents/GBS/Imputation/";
//       dir= "/home/local/MAIZE/kls283/GBS/AllZea2.7MergeWithDepth/";
       dir= "/home/local/MAIZE/kls283/GBS/Imputation/AllZea2.7Genos/";
       String[] files= new String[27];
       int index= 0;
       for (int part = 14; part < 17; part++) {
           if (part<10) files[index]= dir+"part0"+part+"/AllZeaGBS_v2.7_SeqToGenos_part0"+part+".hmp.h5";
           else files[index]= dir+"part"+part+"/AllZeaGBS_v2.7_SeqToGenos_part"+part+".hmp.h5";
           index++;
       }
           ExportUtils.addTaxaFromExistingByteHDF5File(files, dir+"AllZeaGBS_v2.7_SeqToGenos_combined11_14-17.hmp.h5",true);
   }
}
