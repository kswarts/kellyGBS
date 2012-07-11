/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */



import java.io.BufferedOutputStream;
import java.io.BufferedWriter;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.Date;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByteHDF5TagGroups;
import net.maizegenetics.genome.BaseEncoder;
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

       for (int tag= 0; tag < theTBT.getTagCount(); tag++) {//adds the tags that appear in greater than or equal to the taxa specified in the constructor
           if (theTBT.getNumberOfTaxaWithTag(tag) > minSharedTaxa) {
               outputTags.addReadCount(theTBT.getTag(tag), theTBT.getTagLength(tag), 1);
           }
       }
       outputTags.shrinkToCurrentRows();

       TagsByTaxaByte theTBTOut = new TagsByTaxaByte(theTaxa, outputTags);//new blank TBT.Byte file with the taxa from theTBT and the filtered tags
       for (int tagFill=0; tagFill < theTBTOut.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
           for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
               theTBTOut.addReadsToTagTaxon(tagFill, taxonFill, theTBT.getReadCountForTagTaxon(theTBT.getTagIndex(theTBTOut.getTag(tagFill)), taxonFill));
           }
       }
       theTBTOut.writeDistFile(new File(outTBTFileName), FilePacking.Byte, 1);
   }
   
   public static void DissectTOPM(String inTOPMFileName) {
       byte[] byteTest = new byte[1];
       byteTest[0] = Byte.MIN_VALUE;
       System.out.println(byteTest[0]);
       File outTextFileName = new File(inTOPMFileName + ".txt");
       TagsOnPhysicalMap theTOPM = new TagsOnPhysicalMap(inTOPMFileName, true);
       int[] dist = theTOPM.mappingDistribution();

       System.out.println("\nTotal tags:\t" + theTOPM.tagNum
               + "\nMapped tags:\t" + theTOPM.mappedTags()              
               + "\nTag Distribution:\n");
       for (int dis = 0; dis < dist.length; dis++) {
           System.out.println("Tags mapping to " + dis + " locations:\t" + dist[dis]);
       }
       theTOPM.writeTextFile(outTextFileName);


//       int noMap= 0;
//       int oneMap= 0;
//       int twoMap= 0;
//       int greaterTwoMap= 0;
//       a
//       for (int tag= 0; tag < theTOPM.tagNum; tag++) {
//           theTOPM.getPositionArray(tag);
//       }
   }
   /**for biologically characterizing a TOPM and associated tagsByTaxaFile. Positions determined based on AGP_v1, centromeres taken from Wolfgruber et al 2012 and pericentromeric/telomeric extrapolated from Gore et al 2009**/ 
   public static void CharTOPMWithTBT(String inTOPMFileName, String inTBTFileName, int divisor) {
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads in the TOPM, true if binary
       TagsByTaxaByteHDF5TagGroups theTBT= new TagsByTaxaByteHDF5TagGroups(inTBTFileName);
       File outputFile= new File(dir+"/CharTOPMWithTBTLog_unitsOf"+divisor+"TagsFromTBT"+theDateString);
       
       //order TBT by read count
       int[][] orderByRead= new int[theTBT.getTagCount()][2];//a 2D int array to hold tag indices from the topm and the corresponding read counts
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
               taxaNumOfAlignedTags+= theTBT.getNumberOfTaxaWithTag(realTBTIndex); //num taxa that share aligned tag
               stats[division][chr+6]++;//plus 6 because 6 other variables held in array (for chromosomes 0-12)
               stats[division][1]++;//num of tags that align
               
               if (chr > 0 && chr < 11) {
                   if ((minPos < chrPos[1][chr-1]) && (minPos > chrPos[0][chr-1])) stats[division][5]++;                     
                   else if (minPos < chrPos[3][chr-1] && minPos > chrPos[2][chr-1]) stats[division][4]++;                     
                   else stats[division][3]++;
               }
           }
           else {
               taxaNumOfUnalignedTags+= theTBT.getNumberOfTaxaWithTag(realTBTIndex);//num taxa that share unaligned tag
               notAlign++;
           }
           
           if (theTOPM.getPositionArray(TOPMIndex).length > 1) stats[division][2]++;
           
           if (((orderedTBTIndex)%divisor == 0 && orderedTBTIndex != 0) || theTBT.getTagCount()-1 == orderedTBTIndex) {
               taxaStats[division][0]= ((double)taxaNumOfAlignedTags/(double)(align*theTBT.getTaxaCount()))*(double)theTBT.getTaxaCount();//the average number of taxa that share an aligned tag in the divisor set
               taxaStats[division][1]= ((double)taxaNumOfUnalignedTags/(double)(notAlign*theTBT.getTaxaCount()))*(double)theTBT.getTaxaCount();//the average number of taxa that share an unaligned tag in the divisor set
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
       int[][] orderByRead= new int[theTags.getTagCount()][2];//a 2D int array to hold tag indices from the topm and the corresponding read counts
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
   
   public static void main (String args[]) {
       
       dir= "/Users/kelly/Documents/GBS";
//       //args for FilterMergedTBT
//       int minSharedTaxa= 2;
//       String inTBTFileName= dir+ "/landraceCombo/mergedTBT/SW_RI_Span_landrace_min1.tbt.byte";//input tbt file pathway
//       String outTBTFileName= dir+ "/landraceCombo/mergedTBT/SW_RI_Span_landrace_min"+minSharedTaxa+".tbt.byte";//output tbt file
//       FilterMergedTBT(minSharedTaxa, inTBTFileName, outTBTFileName);
       
       //args for CharTOPM
       String inTOPMFileName= dir+"/landraceCombo/topm/SW_RI_Span_landraces_min2_bowtie2MissingOrganelles.topm";
       String inTBTFileName= dir+"/landraceCombo/mergedTBT/SW_RI_Span_landrace_min1.tbt.byte";
       int divisor= 1000000;
       CharTOPMWithTBT(inTOPMFileName, inTBTFileName, divisor);
       
//       //args for CharTOPMWithTagCount
//       String inTOPMFileName= dir+"/AllZeaMasterTags_c10_20120703.topm";
//       String inTagCountFileName= dir+"/AllZeaMasterTags_c10_20120606.cnt";
//       int divisor= 1000000;
//       CharTOPMWithTagCount(inTOPMFileName, inTagCountFileName, divisor);

   }
}