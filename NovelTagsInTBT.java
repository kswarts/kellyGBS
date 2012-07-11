/*
* To change this template, choose Tools | Templates
* and open the template in the editor.
*/



import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCountMutable;
import net.maizegenetics.genome.BaseEncoder;

/** This class contains methods for generating diversity metrics for GBS data, only using tag counts (but standardized by good, barcoded reads)
*
* @author kelly, Jeff Glaubitz
*/
/**This method generates system output for estimating novel genetic diversity by taxa from a tbt file relative to a master tag file, and relative
* to a topm file.  The output types unbiased by B73 are the total number of tags in the input tbt (totalNTags), the number of tags in
* the tbt that do not match a tag present in the input master tag list (nNovel tags), the total number of good barcoded reads present in the tbt
* (totalNGBReads), and the number of good barcoded reads present in the tbt from tags that do not match the input master tag list.  This can be used
* to get the proportion of overall novel genetic variation in your tbt vs the input master tags - using nReadsOfNovelTags/totalNGBReads standardizes
* for different overall read numbers for different taxa. The class also generates output for the number of novel tags that align to a unique location
* in the  B73 reference genome (these have a startPosition > 0) and the number of good barcoded reads that these tags contain from the tbt so that
* the GBR are single copy.  These variables are called nNovelAlignedTagsPerTaxon and nReadsOfNovelAlignedTagsPerTaxon. This should give a better
* estimate of the proportion of novel genic variation found in the samples of the input tbt because most of the tags that don't align are in
* repetitive regions (or that region at the end of chromosome 10(?) that is missing in B73).
*
* Note: The TOPM must contain the same Tags as the TBT!!**/
public class NovelTagsInTBT {
   public static String dir= "//home/local/MAIZE/kls283/GBS";
   public static String compareFileID;
   public static String masterFileID;

   public static void NovelTagsInTBT(){
       String inTBTFileName=        dir+"SW_landraces/mergedTBT/SW_landracesTaxaFilterMin2.tbt.byte";//input tbt file pathway
       String inMasterTagsFileName= dir+"landraceBuild/build20120110MergedTags.cnt";//input a master tag file pathway
       String inTOPMFileName= dir+"SW_landraces/topm/SW_landraces.topm";//input topm file pathway
       TagsByTaxa theTBT= new TagsByTaxaByte(inTBTFileName, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       String[] theTaxa = theTBT.getTaxaNames();//holds taxa names of files in TBT

       int[] totalTagCountPerTaxon = new int[theTaxa.length]; //these arrays hold counts in the same index number as the originating taxa in theTaxa
       int[] nNovelTagsPerTaxon = new int[theTaxa.length];
       int[] totalNGBReadsPerTaxon = new int[theTaxa.length];
       int[] nReadsOfNovelTagsPerTaxon = new int[theTaxa.length];
       int[] nNovelAlignedTagsPerTaxon= new int[theTaxa.length];
       int[] nReadsOfNovelAlignedTagsPerTaxon= new int[theTaxa.length];//single copy tags

       TagCounts myMasterTags = new TagCounts(inMasterTagsFileName, FilePacking.Bit);//reads in the master tag file specified
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads the topm, true if file binary, false if txt. The TOPM should contain the same tags as the TBT

       for (int tag = 0; tag < theTBT.getTagCount(); tag++) { //this loops through the tags in the tbt
           long[] currTag = theTBT.getTag(tag); //This gets the sequence of the current tag.  Tags are stored in binary format as an array of two longs = 128 bits total for 64 bases (A=00,C=01,G=10,T=11)
           int hitMasterTag = myMasterTags.getTagIndex(currTag); //finds index of the tag in the myMasterTags that matches to current tag (negative if tag is novel)
           int B73StartPosOfTag= theTOPM.getStartPosition(theTOPM.getTagIndex(currTag));//the start position of the current tag in the topm. If the tag did not align to B73 it's negative (but the non-aligned tags should still be present in the TOPM)
           for (int taxon = 0; taxon < theTaxa.length; taxon++) { //loop through each taxon in the tbt for the current tag
               totalNGBReadsPerTaxon[taxon] += theTBT.getReadCountForTagTaxon(tag, taxon);//adds GBR from the tbt for the current taxon index of the current tag
               if (theTBT.getReadCountForTagTaxon(tag, taxon) > 0) {//adds to tag count
                   totalTagCountPerTaxon[taxon]++;
               }
               if (hitMasterTag < 0) {//if the current tag matches a tag in the master tag list
                   nReadsOfNovelTagsPerTaxon[taxon] += theTBT.getReadCountForTagTaxon(tag, taxon);//add GBR count for matching tags
                   if (theTBT.getReadCountForTagTaxon(tag, taxon) > 0) nNovelTagsPerTaxon[taxon]++;
                   if (B73StartPosOfTag >= 0) { // if novel tag aligns to unique position in the B73 reference genome
                       if (theTBT.getReadCountForTagTaxon(tag, taxon) >0) {
                           nNovelAlignedTagsPerTaxon[taxon]++;
                           nReadsOfNovelAlignedTagsPerTaxon[taxon]+= theTBT.getReadCountForTagTaxon(tag, taxon);//single copy
                       }
                   }
               }
           }
       }
       System.out.println("Taxon\ttotalNtags\tnNovelTags\tnNovelAlignedTags\ttotalNGBReads\tnGBReadsOfNovelTags\tnGBReadsOfNovelAlignedTags(SingleCopy)");
       for (int taxon = 0; taxon < theTaxa.length; taxon++) {
           System.out.println(theTaxa[taxon]+"\t"+totalTagCountPerTaxon[taxon]+"\t"+nNovelTagsPerTaxon[taxon]+"\t"+nNovelAlignedTagsPerTaxon[taxon]+
                   "\t"+totalNGBReadsPerTaxon[taxon]+"\t"+nReadsOfNovelTagsPerTaxon[taxon]+"\t"+nReadsOfNovelAlignedTagsPerTaxon[taxon]);
       }
   }

   public static void CDFOrderedTagCounts(){
       String inTBTFileNameOne= dir+"/landraceCombo/mergedTBT/SW_RI_Span_landrace_min2.tbt.byte";//use here for landraces
       String inTagCountFileNameTwo= dir+"/B73/mergedTagCounts/B73_min10.cnt";
       String outTBTFileName= dir+"/landraceCombo/mergedTBT/SW_RI_Span_landrace_min2_OrderedByMasterTagReadUbiquity.tbt.byte";
       String inMasterTagsFileName= dir+"/landraceBuild/AllZeaMasterTags_c10_20120606.cnt";//input a master tag file pathway
       String inTOPMFileName= dir+"/landraceBuild/AllZeaMasterTags_c10_20120607.topm";//input topm file pathway
       TagsByTaxa theTBTOne= new TagsByTaxaByte(inTBTFileNameOne, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       TagCounts theTBTTwo= new TagCounts(inTagCountFileNameTwo, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       TagsByTaxa orderedMatchingTBTOne;
       String[] theTaxa = theTBTOne.getTaxaNames();//holds taxa names of files in TBT

       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads the topm, true if file binary, false if txt. The TOPM should contain the same tags as the TBT
       TagCounts myMasterTags = new TagCounts(inMasterTagsFileName, FilePacking.Bit);//reads in the master tag file specified
       TagCountMutable orderMatchingTagOne= new TagCountMutable(theTBTOne.getTagSizeInLong(), theTBTOne.getTagCount());//tag count mutable to hold tags that match TBTOne, ordered by master output
       int[][] orderedMasterTags= new int[myMasterTags.getSize()][2];//a 2D int array to hold tag indices from the master tag list and the corresponding read counts

       for (int masterTag= 0; masterTag < myMasterTags.getSize(); masterTag++) {//fill orderedMasterTags with indices (index 1) and read counts (index 0)
           orderedMasterTags[masterTag][0]= myMasterTags.getReadCount(masterTag);
           orderedMasterTags[masterTag][1]= masterTag;
       }

       Arrays.sort(orderedMasterTags, new java.util.Comparator<int[]>() {//sort 2d array by read count, descending
           @Override
           public int compare(int[] o1, int[] o2) {
               if(o1[0]>o2[0]) return -1; // -1 for descending order
               else if(o1[0]<o2[0]) return 1; // 1 for descending order
               else return 0;
           }
       });

       long[] currTag;
       int hitTagOne;//finds index of the tag in theTBTOne that matches to current tag (negative if tag is not present)
       int hitTagTwo;//finds index of the tag in theTBTTwo that matches to current tag (negative if tag is not present)
       int tagOnePresent= 0;
       int tagTwoPresent= 0;
       int tagOneAlign= 0;
       int tagTwoAlign= 0;
       int masterTagAlign= 0;
       long readCount= 0;
       long alignedReadCount= 0;
       int refStartPosOfTag;
       int millTag= 1;

       System.out.println("Millionth Tag\tMatchingTagOne\tMatchingAlignedTagOne\tMatchingTagTwo\tMatchingAlignedTagTwo\tMasterTagAlignment\tReadCount\tAlignedReadCount");
       for (int masterIndexRead= 0; masterIndexRead < myMasterTags.getSize(); masterIndexRead++) {
           
           if (masterIndexRead%1000000 == 0 && masterIndexRead != 0) {
               System.out.println(millTag+"\t"+tagOnePresent+"\t"+tagOneAlign+"\t"+tagTwoPresent+"\t"+tagTwoAlign+"\t"+masterTagAlign+"\t"+readCount+"\t"+alignedReadCount);
               millTag++;
           }
           
           currTag= myMasterTags.getTag(orderedMasterTags[masterIndexRead][1]);
           hitTagOne= theTBTOne.getTagIndex(currTag);//negative if not present
           hitTagTwo= theTBTTwo.getTagIndex(currTag);//negative if not present
           refStartPosOfTag= theTOPM.getStartPosition(theTOPM.getTagIndex(currTag));//negative if not present
           readCount+= (long) orderedMasterTags[masterIndexRead][0];       
           if (refStartPosOfTag >= 0) {
               masterTagAlign++;
               alignedReadCount+= (long) orderedMasterTags[masterIndexRead][0];
           }
           if (hitTagOne >= 0) {
               tagOnePresent++;
               orderMatchingTagOne.addReadCount(currTag, theTBTOne.getTagLength(theTBTOne.getTagIndex(currTag)), 1);
               if (refStartPosOfTag >= 0) tagOneAlign++;
           }
           if (hitTagTwo >= 0) {
               tagTwoPresent++;
               if (refStartPosOfTag >= 0) tagTwoAlign++;
           }           
       }

       orderMatchingTagOne.collapseCounts();
       orderMatchingTagOne.shrinkToCurrentRows();
       orderedMatchingTBTOne = new TagsByTaxaByte(theTaxa, orderMatchingTagOne);//new blank TBT.Byte file with the taxa from theTBT and the filtered tags
           for (int tagFill=0; tagFill < orderedMatchingTBTOne.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
               for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
                   orderedMatchingTBTOne.addReadsToTagTaxon(tagFill, taxonFill, theTBTOne.getReadCountForTagTaxon(theTBTOne.getTagIndex(orderMatchingTagOne.getTag(tagFill)), taxonFill));
               }
           }
           orderedMatchingTBTOne.writeDistFile(new File(outTBTFileName), FilePacking.Byte, 1);
           KellyUtils.TagsByTaxaToFastq(outTBTFileName, outTBTFileName+".fq");
   }



   public static void NovelTagsInTBTOutputToTBT(int readsPerTagTaxon, int taxaPerTag){
       String fileID= "C08JYACXX_2_min1";
       String inTBTFileName=        dir+"RI_landraces/tbt/"+fileID+".tbt.byte";//input tbt file pathway
       String inMasterTagsFileName= dir+"landraceBuild/build20120110MergedTags.cnt";//input a master tag file pathway
       String inTOPMFileName= dir+"RI_landraces/topm/RI_landraces.topm.bin";//input topm file pathway
       String outAllNovelFileName= dir+"RI_landraces/output/NovelTagsFrom"+fileID+"_minReads"+readsPerTagTaxon+"_minTaxa"+taxaPerTag+".tbt.byte";
       String outUnalignedFileName= dir+"RI_landraces/output/UnalignedNovelTagsFrom"+fileID+"_minReads"+readsPerTagTaxon+"_minTaxa"+taxaPerTag+".tbt.byte";
       String outAlignedFileName= dir+"RI_landraces/output/AlignedNovelTagsFrom"+fileID+"_minReads"+readsPerTagTaxon+"_minTaxa"+taxaPerTag+".tbt.byte";
       File novelTagSummary= new File(dir+"RI_landraces/output/NovelTagsSummaryFrom"+fileID+"_minReads"+readsPerTagTaxon+"_minTaxa"+taxaPerTag+".txt");
       TagsByTaxa theTBT= new TagsByTaxaByte(inTBTFileName, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       String[] theTaxa = theTBT.getTaxaNames();//holds taxa names of files in TBT
       TagsByTaxa allNovelTBTOut;
       TagsByTaxa alignedTBTOut;
       TagsByTaxa unalignedTBTOut;
       TagCountMutable allNovelTags= new TagCountMutable(theTBT.getTagSizeInLong(), theTBT.getTagCount());
       TagCountMutable novelAlignedTags= new TagCountMutable(theTBT.getTagSizeInLong(), theTBT.getTagCount());
       TagCountMutable novelUnalignedTags= new TagCountMutable(theTBT.getTagSizeInLong(), theTBT.getTagCount());

       int[] totalTagCountPerTaxon = new int[theTaxa.length]; //these arrays hold counts in the same index number as the originating taxa in theTaxa
       int[] nNovelTagsPerTaxon = new int[theTaxa.length];
       int[] totalNGBReadsPerTaxon = new int[theTaxa.length];
       int[] nReadsOfNovelTagsPerTaxon = new int[theTaxa.length];
       int[] nNovelAlignedTagsPerTaxon= new int[theTaxa.length];
       int[] nReadsOfNovelAlignedTagsPerTaxon= new int[theTaxa.length];//single copy tags
       int addAligned= 0;
       int addUnaligned= 0;

       TagCounts myMasterTags = new TagCounts(inMasterTagsFileName, FilePacking.Bit);//reads in the master tag file specified
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads the topm, true if file binary, false if txt. The TOPM should contain the same tags as the TBT

       for (int tag = 0; tag < theTBT.getTagCount(); tag++) { //this loops through the tags in the tbt
           long[] currTag = theTBT.getTag(tag); //This gets the sequence of the current tag.  Tags are stored in binary format as an array of two longs = 128 bits total for 64 bases (A=00,C=01,G=10,T=11)
           int hitMasterTag = myMasterTags.getTagIndex(currTag); //finds index of the tag in the myMasterTags that matches to current tag (negative if tag is novel)
           int B73StartPosOfTag= theTOPM.getStartPosition(theTOPM.getTagIndex(currTag));//the start position of the current tag in the topm. If the tag did not align to B73 it's negative (but the non-aligned tags should still be present in the TOPM)
           for (int taxon = 0; taxon < theTaxa.length; taxon++) { //loop through each taxon in the tbt for the current tag
               totalNGBReadsPerTaxon[taxon] += theTBT.getReadCountForTagTaxon(tag, taxon);//adds GBR from the tbt for the current taxon index of the current tag
               if (theTBT.getReadCountForTagTaxon(tag, taxon) > 0) {//adds to tag count
                   totalTagCountPerTaxon[taxon]++;
               }
               if (hitMasterTag < 0) {//if the current tag does not match a tag in the master tag list
                   nReadsOfNovelTagsPerTaxon[taxon] += theTBT.getReadCountForTagTaxon(tag, taxon);//add GBR count for matching tags
                   if (theTBT.getReadCountForTagTaxon(tag, taxon) > 0) nNovelTagsPerTaxon[taxon]++;
                   if (B73StartPosOfTag >= 0) { // if novel tag aligns to unique position in the B73 reference genome
                       if (theTBT.getReadCountForTagTaxon(tag, taxon) >0) {
                           nNovelAlignedTagsPerTaxon[taxon]++;
                           nReadsOfNovelAlignedTagsPerTaxon[taxon]+= theTBT.getReadCountForTagTaxon(tag, taxon);//single copy
                           if (theTBT.getReadCountForTagTaxon(tag, taxon)>readsPerTagTaxon && theTBT.getNumberOfTaxaWithTag(taxon)>taxaPerTag) addAligned++;
                       }
                   }
                   else if (theTBT.getReadCountForTagTaxon(tag, taxon)>readsPerTagTaxon && theTBT.getNumberOfTaxaWithTag(taxon)>taxaPerTag) addUnaligned++;
               }
           }
           if(addUnaligned>0) novelUnalignedTags.addReadCount(currTag, theTBT.getTagLength(tag), 1);
           if(addAligned>0) novelAlignedTags.addReadCount(currTag, theTBT.getTagLength(tag), 1);
           if(addUnaligned>0 ||addAligned>0) allNovelTags.addReadCount(currTag, theTBT.getTagLength(tag), 1);
           addUnaligned= 0;
           addAligned= 0;
       }

       System.out.println(novelAlignedTags.getCurrentSize());
       allNovelTags.shrinkToCurrentRows();
       novelAlignedTags.shrinkToCurrentRows();
       novelUnalignedTags.shrinkToCurrentRows();
       System.out.println(novelAlignedTags.getCurrentSize());

       System.out.println("Taxon\ttotalNtags\tnNovelTags\tnNovelAlignedTags\ttotalNGBReads\tnGBReadsOfNovelTags\tnGBReadsOfNovelAlignedTags(SingleCopy)");
       for (int taxon = 0; taxon < theTaxa.length; taxon++) {
           System.out.println(theTaxa[taxon]+"\t"+totalTagCountPerTaxon[taxon]+"\t"+nNovelTagsPerTaxon[taxon]+"\t"+nNovelAlignedTagsPerTaxon[taxon]+
                   "\t"+totalNGBReadsPerTaxon[taxon]+"\t"+nReadsOfNovelTagsPerTaxon[taxon]+"\t"+nReadsOfNovelAlignedTagsPerTaxon[taxon]);
       }

       allNovelTBTOut = new TagsByTaxaByte(theTaxa, allNovelTags);//new blank TBT.Byte file with the taxa from theTBT and the filtered tags
       for (int tagFill=0; tagFill < allNovelTBTOut.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
           for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
               allNovelTBTOut.addReadsToTagTaxon(tagFill, taxonFill, theTBT.getReadCountForTagTaxon(theTBT.getTagIndex(allNovelTags.getTag(tagFill)), taxonFill));
           }
       }
       allNovelTBTOut.writeDistFile(new File(outAllNovelFileName), FilePacking.Byte, 1);

       unalignedTBTOut = new TagsByTaxaByte(theTaxa, novelUnalignedTags);//new blank TBT.Byte file with the taxa from theTBT and the filtered tags
       for (int tagFill=0; tagFill < unalignedTBTOut.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
           for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
               unalignedTBTOut.addReadsToTagTaxon(tagFill, taxonFill, theTBT.getReadCountForTagTaxon(theTBT.getTagIndex(novelUnalignedTags.getTag(tagFill)), taxonFill));
           }
       }
       unalignedTBTOut.writeDistFile(new File(outUnalignedFileName), FilePacking.Byte, 1);

       alignedTBTOut = new TagsByTaxaByte(theTaxa, novelAlignedTags);//new blank TBT.Byte file with the taxa from theTBT and the filtered tags
       for (int tagFill=0; tagFill < novelAlignedTags.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
           for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
               alignedTBTOut.addReadsToTagTaxon(tagFill, taxonFill, theTBT.getReadCountForTagTaxon(theTBT.getTagIndex(novelAlignedTags.getTag(tagFill)), taxonFill));
           }
       }
       alignedTBTOut.writeDistFile(new File(outAlignedFileName), FilePacking.Byte, 1);

       try {
           novelTagSummary.createNewFile();
           novelTagSummary.setWritable(true);
           BufferedWriter BW = new BufferedWriter(new FileWriter(novelTagSummary));

           BW.write("TagSequence\tchr\tstrand\tstartPos\n");
           int[] posArray;
           for (int h=0; h < allNovelTags.getTotalCount(); h++) {
               BW.write(BaseEncoder.getSequenceFromLong(allNovelTags.getTag(h)));
               if (theTOPM.getStartPosition(theTOPM.getTagIndex(allNovelTags.getTag(h)))>0) {
                   posArray= theTOPM.getPositionArray(theTOPM.getTagIndex(allNovelTags.getTag(h)));
                   BW.write("\t"+posArray[0]+"\t"+posArray[1]+"\t"+posArray[2]+"\n");
               }
               else BW.write("\n");
           }
           BW.close();
       }

       catch (IOException e) {
           System.out.println(e);
       }
   }

   public static void NovelTagsInTBTOutputToTBTReadsPerTag(int[] readsPerTag, int[] taxaPerTag){
       String localDirOut= "/landraceCombo/output/NovelTagsInTBTOutputToTBTAgainstJulyBuildBowtie2/";
       String localDirIn= "/landraceCombo/";
       String fileID= "/SW_RI_Span_landrace_min1";
       String inTBTFileName=        dir+localDirIn+"mergedTBT/"+fileID+".tbt.byte";//input tbt file pathway
//        String inMasterTagsFileName= dir+"landraceBuild/build20120110MergedTags.cnt";//input a master tag file pathway
//        String inTOPMFileName= dir+localDirIn+"topm/SW_RI_Span_landraces_min2.topm";//input topm file pathway
//        String fileID= "SW_RI_Span_landrace_min1";
//        String inTBTFileName=        dir+localDirIn+"mergedTBT/"+fileID+".tbt.byte";//input tbt file pathway
       String inMasterTagsFileName= dir+"/landraceBuild/AllZeaMasterTags_c10_20120606.cnt";//input a master tag file pathway
       String inTOPMFileName= dir+"/landraceCombo/topm/AllZeaMasterTags_c10_20120607.topm";//input topm file pathway
//        String fileID= "SW_landraces_min1";
//        String inTBTFileName=        dir+"SW_landraces/mergedTBT/"+fileID+".tbt.byte";//input tbt file pathway
//        String inMasterTagsFileName= dir+"landraceBuild/build20120110MergedTags.cnt";//input a master tag file pathway
//        String inTOPMFileName= dir+"SW_landraces/topm/SW_landracesMinRead2.topm";//input topm file pathway

       TagsByTaxa theTBT= new TagsByTaxaByte(inTBTFileName, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       TagsByTaxa novelTags;
       TagCounts theTags= new TagCounts(inTBTFileName, FilePacking.Byte);
       String[] theTaxa = theTBT.getTaxaNames();//holds taxa names of files in TBT
       TagsByTaxa allNovelTBTOut;
       TagsByTaxa alignedTBTOut;
       TagsByTaxa unalignedTBTOut;
       TagCounts myMasterTags = new TagCounts(inMasterTagsFileName, FilePacking.Bit);//reads in the master tag file specified
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads the topm, true if file binary, false if txt. The TOPM should contain the same tags as the TBT

       for (int minTag= 0; minTag < readsPerTag.length; minTag++) {
           for (int minTaxa= 0; minTaxa < taxaPerTag.length; minTaxa++) {
               System.out.println("\nBegin new run with filters min reads per tag: "+readsPerTag[minTag]+" and min taxa per tag: "+taxaPerTag[minTaxa]);
               String outAllNovelFileName= dir+localDirOut+"NovelTagsReadPerTagFrom"+fileID+"_minReads"+readsPerTag[minTag]+"_minTaxa"+taxaPerTag[minTaxa]+".tbt.byte";
               String outUnalignedFileName= dir+localDirOut+"UnalignedNovelTagsReadPerTagFrom"+fileID+"_minReads"+readsPerTag[minTag]+"_minTaxa"+taxaPerTag[minTaxa]+".tbt.byte";
               String outAlignedFileName= dir+localDirOut+"AlignedNovelTagsReadPerTagFrom"+fileID+"_minReads"+readsPerTag[minTag]+"_minTaxa"+taxaPerTag[minTaxa]+".tbt.byte";
               File novelTagSummary= new File(dir+localDirOut+"NovelTagsReadPerTagSummaryFrom"+fileID+"_minReads"+readsPerTag[minTag]+"_minTaxa"+taxaPerTag[minTaxa]+".txt");
               TagCountMutable allNovelTags= new TagCountMutable(theTBT.getTagSizeInLong(), theTBT.getTagCount());
               TagCountMutable novelAlignedTags= new TagCountMutable(theTBT.getTagSizeInLong(), theTBT.getTagCount());
               TagCountMutable novelUnalignedTags= new TagCountMutable(theTBT.getTagSizeInLong(), theTBT.getTagCount());
               int[] totalTagCountPerTaxon = new int[theTaxa.length]; //these arrays hold counts in the same index number as the originating taxa in theTaxa
               int[] nNovelTagsPerTaxon = new int[theTaxa.length];
               int[] totalNGBReadsPerTaxon = new int[theTaxa.length];
               int[] nReadsOfNovelTagsPerTaxon = new int[theTaxa.length];
               int[] nNovelAlignedTagsPerTaxon= new int[theTaxa.length];
               int[] nReadsOfNovelAlignedTagsPerTaxon= new int[theTaxa.length];//single copy tags
               int addAligned= 0;
               int addUnaligned= 0;

               for (int tag = 0; tag < theTBT.getTagCount(); tag++) { //this loops through the tags in the tbt
                   long[] currTag = theTBT.getTag(tag); //This gets the sequence of the current tag.  Tags are stored in binary format as an array of two longs = 128 bits total for 64 bases (A=00,C=01,G=10,T=11)
                   int hitMasterTag = myMasterTags.getTagIndex(currTag); //finds index of the tag in the myMasterTags that matches to current tag (negative if tag is novel)
                   int B73StartPosOfTag= theTOPM.getStartPosition(theTOPM.getTagIndex(currTag));//the start position of the current tag in the topm. If the tag did not align to B73 it's negative (but the non-aligned tags should still be present in the TOPM)
                   for (int taxon = 0; taxon < theTaxa.length; taxon++) { //loop through each taxon in the tbt for the current tag
                       totalNGBReadsPerTaxon[taxon] += theTBT.getReadCountForTagTaxon(tag, taxon);//adds GBR from the tbt for the current taxon index of the current tag
                       if (theTBT.getReadCountForTagTaxon(tag, taxon) > 0) {//adds to tag count
                           totalTagCountPerTaxon[taxon]++;
                       }
                       if (hitMasterTag < 0) {//if the current tag does not match a tag in the master tag list
                           nReadsOfNovelTagsPerTaxon[taxon] += theTBT.getReadCountForTagTaxon(tag, taxon);//add GBR count for matching tags
                           if (theTBT.getReadCountForTagTaxon(tag, taxon) > 0) nNovelTagsPerTaxon[taxon]++;
                           if (B73StartPosOfTag >= 0) { // if novel tag aligns to unique position in the B73 reference genome
                               if (theTBT.getReadCountForTagTaxon(tag, taxon) >0) {
                                   nNovelAlignedTagsPerTaxon[taxon]++;
                                   nReadsOfNovelAlignedTagsPerTaxon[taxon]+= theTBT.getReadCountForTagTaxon(tag, taxon);//single copy
                                   if (theTags.getReadCount(tag)>readsPerTag[minTag] && theTBT.getNumberOfTaxaWithTag(taxon)>taxaPerTag[minTaxa]) addAligned++;
                               }
                           }
                           else if (theTags.getReadCount(tag)>readsPerTag[minTag] && theTBT.getNumberOfTaxaWithTag(taxon)>taxaPerTag[minTaxa]) addUnaligned++;
                       }
                   }
                   if(addUnaligned>0) novelUnalignedTags.addReadCount(currTag, theTBT.getTagLength(tag), 1);
                   if(addAligned>0) novelAlignedTags.addReadCount(currTag, theTBT.getTagLength(tag), 1);
                   if(addUnaligned>0 ||addAligned>0) allNovelTags.addReadCount(currTag, theTBT.getTagLength(tag), 1);
                   addUnaligned= 0;
                   addAligned= 0;
               }

               System.out.println(novelAlignedTags.getCurrentSize());
               allNovelTags.shrinkToCurrentRows();
               novelAlignedTags.shrinkToCurrentRows();
               novelUnalignedTags.shrinkToCurrentRows();
               System.out.println(novelAlignedTags.getCurrentSize());

               System.out.println("Taxon\ttotalNtags\tnNovelTags\tnNovelAlignedTags\ttotalNGBReads\tnGBReadsOfNovelTags\tnGBReadsOfNovelAlignedTags(SingleCopy)");
               for (int taxon = 0; taxon < theTaxa.length; taxon++) {
                   System.out.println(theTaxa[taxon]+"\t"+totalTagCountPerTaxon[taxon]+"\t"+nNovelTagsPerTaxon[taxon]+"\t"+nNovelAlignedTagsPerTaxon[taxon]+
                           "\t"+totalNGBReadsPerTaxon[taxon]+"\t"+nReadsOfNovelTagsPerTaxon[taxon]+"\t"+nReadsOfNovelAlignedTagsPerTaxon[taxon]);
               }

               allNovelTBTOut = new TagsByTaxaByte(theTaxa, allNovelTags);//new blank TBT.Byte file with the taxa from theTBT and the filtered tags
               for (int tagFill=0; tagFill < allNovelTBTOut.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
                   for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
                       allNovelTBTOut.addReadsToTagTaxon(tagFill, taxonFill, theTBT.getReadCountForTagTaxon(theTBT.getTagIndex(allNovelTags.getTag(tagFill)), taxonFill));
                   }
               }
               allNovelTBTOut.writeDistFile(new File(outAllNovelFileName), FilePacking.Byte, 1);
               KellyUtils.TagsByTaxaToFastq(outAllNovelFileName, outAllNovelFileName+".fq");

               unalignedTBTOut = new TagsByTaxaByte(theTaxa, novelUnalignedTags);//new blank TBT.Byte file with the taxa from theTBT and the filtered tags
               for (int tagFill=0; tagFill < unalignedTBTOut.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
                   for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
                       unalignedTBTOut.addReadsToTagTaxon(tagFill, taxonFill, theTBT.getReadCountForTagTaxon(theTBT.getTagIndex(novelUnalignedTags.getTag(tagFill)), taxonFill));
                   }
               }
               unalignedTBTOut.writeDistFile(new File(outUnalignedFileName), FilePacking.Byte, 1);
               KellyUtils.TagsByTaxaToFastq(outUnalignedFileName, outUnalignedFileName+".fq");

               alignedTBTOut = new TagsByTaxaByte(theTaxa, novelAlignedTags);//new blank TBT.Byte file with the taxa from theTBT and the filtered tags
               for (int tagFill=0; tagFill < novelAlignedTags.getTagCount(); tagFill++) { //fills the new TBT.Byte with information from theTBT
                   for (int taxonFill=0; taxonFill < theTaxa.length; taxonFill++) {
                       alignedTBTOut.addReadsToTagTaxon(tagFill, taxonFill, theTBT.getReadCountForTagTaxon(theTBT.getTagIndex(novelAlignedTags.getTag(tagFill)), taxonFill));
                   }
               }
               alignedTBTOut.writeDistFile(new File(outAlignedFileName), FilePacking.Byte, 1);
               KellyUtils.TagsByTaxaToFastq(outAlignedFileName, outAlignedFileName+".fq");

//                try {
//                    novelTagSummary.createNewFile();
//                    novelTagSummary.setWritable(true);
//                    BufferedWriter BW = new BufferedWriter(new FileWriter(novelTagSummary));
//
//                    BW.write("TagSequence\tchr\tstrand\tstartPos\n");
//                    int[] posArray;
//                    for (int h=0; h < allNovelTags.getTotalCount(); h++) {
//                        BW.write(BaseEncoder.getSequenceFromLong(allNovelTags.getTag(h)));
//                        if (theTOPM.getPositionMin(theTOPM.getTagIndex(allNovelTags.getTag(h)))>0) {
//                            posArray= theTOPM.getPositionArray(theTOPM.getTagIndex(allNovelTags.getTag(h)));
//                            BW.write("\t"+posArray[0]+"\t"+posArray[1]+"\t"+posArray[2]+"\n");
//                        }
//                        else BW.write("\n");
//                    }
//                    BW.close();
//                }
//
//                catch (IOException e) {
//                    System.out.println(e);
//                }
           }
       }
   }

   /**for each taxon, get tag differences for each taxon relative to each other taxon in TBT**/
   public static void PairwiseTagComparisonsInTBT() {
       File alignedPairwiseFile= new File("/Users/kelly/Documents/GBS/SW_landraces/output/SW_landracesAlignedPairwise_noTagMin.txt");
       File pairwiseFile= new File("/Users/kelly/Documents/GBS/SW_landraces/output/SW_landracespairwise_noTagMin.txt");

       String inTBTFileName= "/Users/kelly/Documents/GBS/SW_landraces/tbt/C0F9UACXX_1_noTagCountMin.tbt.byte";//input tbt file pathway
       String inTOPMFileName= "/Users/kelly/Documents/GBS/SW_landraces/topm/SW_landracesMin1.topm";//input topm file pathway
       TagsByTaxa theTBT= new TagsByTaxaByte(inTBTFileName, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads the topm, true if file binary, false if txt. The TOPM should contain the same tags as the TBT
       String[] theTaxa = theTBT.getTaxaNames();//holds taxa names of files in TBT
       double[][] pairwiseMatrix= new double[theTaxa.length][theTaxa.length]; //holds NGV of taxa in columns vs taxa in rows
       double[][] alignedPairwiseMatrix= new double[theTaxa.length][theTaxa.length]; //holds aligned NGV of taxa in columns vs taxa in rows

       for (int taxa = 0; taxa < theTaxa.length; taxa++) { //loops through the taxa

           for (int taxonMatch = 0; taxonMatch < theTaxa.length; taxonMatch++) { //loops through the taxa to match taxon in pairwise fashion
               int readsInMatchedTagsInCurrTaxon= 0;
               int totalReadsInCurrTaxon= 0;
               int alignedReadsInMatchedTagsInCurrTaxon= 0;
               int alignedTotalReadsInCurrTaxon= 0;

               for (int tag = 0; tag < theTBT.getTagCount(); tag++) { //loops through the tags and adds reads from tags in current taxon that are shared by the comparison taxon and total reads from the current taxon
                   int currTaxonTagCount= theTBT.getReadCountForTagTaxon(tag, taxa);
                   int matchTaxonTagCount= theTBT.getReadCountForTagTaxon(tag, taxonMatch);
                   int B73StartPosOfTag= theTOPM.getStartPosition(theTOPM.getTagIndex(theTBT.getTag(tag)));//the start position of the current tag in the topm. If the tag did not align to B73 it's negative (but the non-aligned tags should still be present in the TOPM)
                   totalReadsInCurrTaxon+= currTaxonTagCount;
                   if (B73StartPosOfTag >= 0) {
                           alignedTotalReadsInCurrTaxon+= currTaxonTagCount;
                       }
                   if (currTaxonTagCount > 0 & matchTaxonTagCount > 0) {
                       readsInMatchedTagsInCurrTaxon+= currTaxonTagCount;
                       if (B73StartPosOfTag >= 0) {
                           alignedReadsInMatchedTagsInCurrTaxon+= currTaxonTagCount;
                       }
                   }
               }
               pairwiseMatrix[taxa][taxonMatch]= (double)readsInMatchedTagsInCurrTaxon/(double)totalReadsInCurrTaxon;
               alignedPairwiseMatrix[taxa][taxonMatch]= (double) alignedReadsInMatchedTagsInCurrTaxon/(double)alignedTotalReadsInCurrTaxon;
//                System.out.println(pairwiseMatrix[taxa][taxonMatch]);
//                System.out.println(alignedPairwiseMatrix[taxa][taxonMatch]);
           }
           System.out.println("finish taxa: " +taxa+"\t"+"taxon name: "+theTBT.getTaxaName(taxa));
       }

       try {

           pairwiseFile.createNewFile();
           pairwiseFile.setWritable(true);
           BufferedWriter pairwiseBW = new BufferedWriter(new FileWriter(pairwiseFile));
           double num= 0.0;

           for (int taxa = 0; taxa < theTaxa.length; taxa++) {

                pairwiseBW.write(theTaxa[taxa]+"\t");
            }
           pairwiseBW.write("\n");

           for (int j = 0; j < theTaxa.length; j++) {
               pairwiseBW.write(theTaxa[j]+"\t");

               for (int i = 0; i < theTaxa.length; i++) {
                   num= pairwiseMatrix[j][i];
                   pairwiseBW.write(Double.toString(num)+"\t");
               }
               pairwiseBW.write("\n");

           }

       pairwiseBW.close();
       }

       catch (IOException e) {
               System.out.println(e);
           }

       try {

            alignedPairwiseFile.createNewFile();
            alignedPairwiseFile.setWritable(true);
            BufferedWriter alignedPairwiseBW = new BufferedWriter(new FileWriter(alignedPairwiseFile));
            double numAligned = 0.0;

            for (int taxa = 0; taxa < theTaxa.length; taxa++) {

                alignedPairwiseBW.write(theTaxa[taxa]+"\t");
            }

            alignedPairwiseBW.write("\n");

            for (int j = 0; j < theTaxa.length; j++) {
                alignedPairwiseBW.write(theTaxa[j]+"\t");

               for (int i = 0; i < theTaxa.length; i++) {
                   numAligned= alignedPairwiseMatrix[j][i];
                   alignedPairwiseBW.write(Double.toString(numAligned)+"\t");
               }
               alignedPairwiseBW.write("\n");

           }

       alignedPairwiseBW.close();
       }

       catch (IOException e) {
               System.out.println(e);
           }
   }

   public static void PairwiseTagComparisonsWithSharedReadCorrection() {
       File alignedPairwiseFile= new File("/Users/kelly/Documents/GBS/SW_landraces/output/SW_landracesAlignedPairwise_readCorrected.txt");
       File pairwiseFile= new File("/Users/kelly/Documents/GBS/SW_landraces/output/SW_landracespairwise_readCorrected.txt");

       String inTBTFileName= "/Users/kelly/Documents/GBS/SW_landraces/mergedTBT/SW_landraces_wSorghum_TaxaFilterMin2.tbt.byte";//input tbt file pathway
       String inTOPMFileName= "/Users/kelly/Documents/GBS/SW_landraces/topm/SW_landraces_wSB.topm";//input topm file pathway
       TagsByTaxa theTBT= new TagsByTaxaByte(inTBTFileName, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads the topm, true if file binary, false if txt. The TOPM should contain the same tags as the TBT
       String[] theTaxa= theTBT.getTaxaNames();//holds taxa names of files in TBT
       int[] theTaxaGoodReadCount= new int[theTaxa.length];//holds the total good read count for each taxon, in the same order as theTaxa
       int[] theTaxaTagCount= new int[theTaxa.length];//holds the total tag count for each taxon, in the same order as theTaxa
       double[][] pairwiseMatrix= new double[theTaxa.length][theTaxa.length]; //holds NGV of taxa in columns vs taxa in rows
       double[][] alignedPairwiseMatrix= new double[theTaxa.length][theTaxa.length]; //holds aligned NGV of taxa in columns vs taxa in rows
       List<Integer> shuffledTagIndices= new ArrayList<Integer>();//holds a list of int

       for(int tagNum = 0; tagNum < theTBT.getTagCount(); tagNum++) {//fills the shuffleTags list with sucessive numbers, up to the number of tags in the TBT
           shuffledTagIndices.add(tagNum);
       }
//        System.out.println("shuffledTagIndices: "+shuffledTagIndices.get(0)+shuffledTagIndices.get(1)+shuffledTagIndices.get(2)+shuffledTagIndices.get(3));

       for(int taxonRead = 0; taxonRead < theTaxa.length; taxonRead++) {//fills theTaxaGoodReadCount and theTaxaTagCount

           for (int tagCount = 0; tagCount < theTBT.getTagCount(); tagCount++) {
               theTaxaGoodReadCount[taxonRead]+= theTBT.getReadCountForTagTaxon(tagCount, taxonRead);
               if (theTBT.getReadCountForTagTaxon(tagCount, taxonRead) > 0) {
                   theTaxaTagCount[taxonRead]++;
               }
           }
//            System.out.println("taxon: "+theTaxa[taxonRead]+"\t"+"taxonReadCount: "+theTaxaGoodReadCount[taxonRead]+"\t"+"taxonTagCount: "+theTaxaTagCount[taxonRead]);
       }
       System.out.println("fill theTaxaGoodReadCount and theTaxaTagCount complete");

       for (int taxa = 0; taxa < theTaxa.length; taxa++) { //loops through the taxa

           for (int taxonMatch = 0; taxonMatch < theTaxa.length; taxonMatch++) { //loops through the taxa to match taxon in pairwise fashion
               int readsInMatchedTagsInCurrTaxon= 0;
               int totalReadsofSampledTagsInCurrTaxon= 0;
               int alignedReadsInMatchedTagsInCurrTaxon= 0;
               int alignedTotalReadsInMatchedTagsInCurrTaxon= 0;
               Collections.shuffle(shuffledTagIndices);
               int tagSampleLimit= (Math.min(theTaxaGoodReadCount[taxa],theTaxaGoodReadCount[taxonMatch])==theTaxaGoodReadCount[taxa]?theTaxaTagCount[taxa]:theTaxaTagCount[taxonMatch]);//figures out which of the two taxa being compared has the lower GBRead and holds the tag count for the lower taxon
//                System.out.println("shuffledTagIndices: "+shuffledTagIndices.get(0)+"\t"+shuffledTagIndices.get(1)+"\t"+shuffledTagIndices.get(2)+"\t"+shuffledTagIndices.get(3));
//                System.out.println("tag sample limit: "+tagSampleLimit);

               for (int tag = 0; tag < tagSampleLimit; tag++) { //samples tags up through the tag count of the taxon with fewer good barcoded reads
                   int currTagIndex= shuffledTagIndices.get(tag);
                   int currTaxonTagCount= theTBT.getReadCountForTagTaxon(currTagIndex, taxa);
                   int matchTaxonTagCount= theTBT.getReadCountForTagTaxon(currTagIndex, taxonMatch);
                   int B73StartPosOfTag= theTOPM.getStartPosition(theTOPM.getTagIndex(theTBT.getTag(currTagIndex)));//the start position of the current tag in the topm. If the tag did not align to B73 it's negative (but the non-aligned tags should still be present in the TOPM)
                   totalReadsofSampledTagsInCurrTaxon+= currTaxonTagCount;
                   if (B73StartPosOfTag >= 0) {
                           alignedTotalReadsInMatchedTagsInCurrTaxon+= currTaxonTagCount;
                       }
                   if (currTaxonTagCount > 0 & matchTaxonTagCount > 0) {
                       readsInMatchedTagsInCurrTaxon+= currTaxonTagCount;
                       if (B73StartPosOfTag >= 0) {
                           alignedReadsInMatchedTagsInCurrTaxon+= currTaxonTagCount;
                       }
                   }
               }
               pairwiseMatrix[taxa][taxonMatch]= (double)readsInMatchedTagsInCurrTaxon/(double)totalReadsofSampledTagsInCurrTaxon;
               alignedPairwiseMatrix[taxa][taxonMatch]= (double) alignedReadsInMatchedTagsInCurrTaxon/(double)alignedTotalReadsInMatchedTagsInCurrTaxon;
//                System.out.println(pairwiseMatrix[taxa][taxonMatch]);
//                System.out.println(alignedPairwiseMatrix[taxa][taxonMatch]);
           }
           System.out.println("finish taxa: " +taxa+"\t"+"taxon name: "+theTBT.getTaxaName(taxa)+"\t"+"taxon good barcoded reads: "+theTaxaGoodReadCount[taxa]+"\t"+"taxon tag count: "+theTaxaTagCount[taxa]);
       }

       try {

           pairwiseFile.createNewFile();
           pairwiseFile.setWritable(true);
           BufferedWriter pairwiseBW = new BufferedWriter(new FileWriter(pairwiseFile));
           double num= 0.0;

           pairwiseBW.write("\t");
           for (int taxa = 0; taxa < theTaxa.length; taxa++) {

                pairwiseBW.write(theTaxa[taxa]+"\t");
            }
           pairwiseBW.write("\n");

           for (int j = 0; j < theTaxa.length; j++) {
               pairwiseBW.write(theTaxa[j]+"\t");

               for (int i = 0; i < theTaxa.length; i++) {
                   num= pairwiseMatrix[j][i];
                   pairwiseBW.write(Double.toString(num)+"\t");
               }
               pairwiseBW.write("\n");

           }

       pairwiseBW.close();
       }

       catch (IOException e) {
               System.out.println(e);
           }

       try {

            alignedPairwiseFile.createNewFile();
            alignedPairwiseFile.setWritable(true);
            BufferedWriter alignedPairwiseBW = new BufferedWriter(new FileWriter(alignedPairwiseFile));
            double numAligned = 0.0;

            alignedPairwiseBW.write("\t");
            for (int taxa = 0; taxa < theTaxa.length; taxa++) {

                alignedPairwiseBW.write(theTaxa[taxa]+"\t");
            }

            alignedPairwiseBW.write("\n");

            for (int j = 0; j < theTaxa.length; j++) {
                alignedPairwiseBW.write(theTaxa[j]+"\t");

               for (int i = 0; i < theTaxa.length; i++) {
                   numAligned= alignedPairwiseMatrix[j][i];
                   alignedPairwiseBW.write(Double.toString(numAligned)+"\t");
               }
               alignedPairwiseBW.write("\n");

           }

       alignedPairwiseBW.close();
       }

       catch (IOException e) {
               System.out.println(e);
           }
   }

   public static void SymmetricalPairwiseTagComparisonsWithSharedReadCorrection() {
       File alignedPairwiseFile= new File("/Users/kelly/Documents/GBS/SW_landraces/output/SW_landracesAlignedPairwise_SymReadCorr.txt");
       File pairwiseFile= new File("/Users/kelly/Documents/GBS/SW_landraces/output/SW_landracespairwise_SymReadCorr.txt");

       String inTBTFileName= "/Users/kelly/Documents/GBS/SW_landraces/mergedTBT/SW_landraces_wSorghum_TaxaFilterMin3.tbt.byte";//input tbt file pathway
       String inTOPMFileName= "/Users/kelly/Documents/GBS/SW_landraces/topm/SW_landraces_wSB.topm";//input topm file pathway
       TagsByTaxa theTBT= new TagsByTaxaByte(inTBTFileName, FilePacking.Byte);//reads in the tbt file specified **modify based on primitive data type
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);//reads the topm, true if file binary, false if txt. The TOPM should contain the same tags as the TBT
       String[] theTaxa= theTBT.getTaxaNames();//holds taxa names of files in TBT
       int[] theTaxaGoodReadCount= new int[theTaxa.length];//holds the total good read count for each taxon, in the same order as theTaxa
       int[] theTaxaTagCount= new int[theTaxa.length];//holds the total tag count for each taxon, in the same order as theTaxa
       int denom= 0;//tracks the GBReads used in the denominator
       double[][] pairwiseMatrix= new double[theTaxa.length][theTaxa.length]; //holds NGV of taxa in columns vs taxa in rows
       double[][] alignedPairwiseMatrix= new double[theTaxa.length][theTaxa.length]; //holds aligned NGV of taxa in columns vs taxa in rows
       List<Integer> ShuffleTagIndices= new ArrayList<Integer>();//holds a list of int

       for(int tagNum = 0; tagNum < theTBT.getTagCount(); tagNum++) {//fills the shuffleTags list with sucessive numbers, up to the number of tags in the TBT
           ShuffleTagIndices.add(tagNum);
       }
//        System.out.println("shuffledTagIndices: "+shuffledTagIndices.get(0)+shuffledTagIndices.get(1)+shuffledTagIndices.get(2)+shuffledTagIndices.get(3));

       for(int taxonRead = 0; taxonRead < theTaxa.length; taxonRead++) {//fills theTaxaGoodReadCount and theTaxaTagCount

           for (int tagCount = 0; tagCount < theTBT.getTagCount(); tagCount++) {
               theTaxaGoodReadCount[taxonRead]+= theTBT.getReadCountForTagTaxon(tagCount, taxonRead);
               if (theTBT.getReadCountForTagTaxon(tagCount, taxonRead) > 0) {
                   theTaxaTagCount[taxonRead]++;
               }
           }
//            System.out.println("taxon: "+theTaxa[taxonRead]+"\t"+"taxonReadCount: "+theTaxaGoodReadCount[taxonRead]+"\t"+"taxonTagCount: "+theTaxaTagCount[taxonRead]);
       }
       System.out.println("fill theTaxaGoodReadCount and theTaxaTagCount complete");

       for (int taxa = 0; taxa < theTaxa.length; taxa++) { //loops through the taxa

           for (int taxonMatch = 0; taxonMatch < theTaxa.length; taxonMatch++) { //loops through the taxa to match taxon in pairwise fashion
               int minPairReadsInMatchedTags= 0;
               int totalReadsofSampledTagsInMinPair= 0;
               int alignedMinPairReadsInMatchedTags= 0;
               int alignedTotalReadsInMatchedTagsInMinPair= 0;
               int pairMin= (Math.min(theTaxaGoodReadCount[taxa],theTaxaGoodReadCount[taxonMatch])==theTaxaGoodReadCount[taxa]?taxa:taxonMatch);//figures out which of the two taxa being compared has the lower GBRead and holds the index value of the lower read in pairMin
               int tagSampleLimit= theTaxaTagCount[pairMin];// holds the tag count for the taxon with fewer reads
               Collections.shuffle(ShuffleTagIndices);
//                System.out.println("shuffledTagIndices: "+shuffledTagIndices.get(0)+"\t"+shuffledTagIndices.get(1)+"\t"+shuffledTagIndices.get(2)+"\t"+shuffledTagIndices.get(3));
//                System.out.println("tag sample limit: "+tagSampleLimit);

               for (int tag = 0; tag < tagSampleLimit; tag++) { //samples tags up through the tag count of the taxon with fewer good barcoded reads
                   int currTagIndex= ShuffleTagIndices.get(tag);
                   int currTaxonTagCount= theTBT.getReadCountForTagTaxon(currTagIndex, taxa);
                   int matchTaxonTagCount= theTBT.getReadCountForTagTaxon(currTagIndex, taxonMatch);
                   int pairMinTagCount= theTBT.getReadCountForTagTaxon(currTagIndex, pairMin);
                   int B73StartPosOfTag= theTOPM.getStartPosition(theTOPM.getTagIndex(theTBT.getTag(currTagIndex)));//the start position of the current tag in the topm. If the tag did not align to B73 it's negative (but the non-aligned tags should still be present in the TOPM)
                   totalReadsofSampledTagsInMinPair+= pairMinTagCount;
                   if (B73StartPosOfTag >= 0) {
                           alignedTotalReadsInMatchedTagsInMinPair+= pairMinTagCount;
                       }
                   if (currTaxonTagCount > 0 & matchTaxonTagCount > 0) {
                       minPairReadsInMatchedTags+= pairMinTagCount;
                       if (B73StartPosOfTag >= 0) {
                           alignedMinPairReadsInMatchedTags+= pairMinTagCount;
                       }
                   }
               }
               pairwiseMatrix[taxa][taxonMatch]= (double)minPairReadsInMatchedTags/(double)totalReadsofSampledTagsInMinPair;
               denom= totalReadsofSampledTagsInMinPair;
               alignedPairwiseMatrix[taxa][taxonMatch]= (double) alignedMinPairReadsInMatchedTags/(double)alignedTotalReadsInMatchedTagsInMinPair;
//                System.out.println(pairwiseMatrix[taxa][taxonMatch]);
//                System.out.println(alignedPairwiseMatrix[taxa][taxonMatch]);
           }
           System.out.println("finish taxa: " +taxa+"\t"+"taxon name: "+theTBT.getTaxaName(taxa)+"\t"+"taxon good barcoded reads: "+theTaxaGoodReadCount[taxa]+"\t"+"taxon tag count: "+theTaxaTagCount[taxa]+"\t"+"denominator: "+denom);
       }

       try {

           pairwiseFile.createNewFile();
           pairwiseFile.setWritable(true);
           BufferedWriter pairwiseBW = new BufferedWriter(new FileWriter(pairwiseFile));
           double num= 0.0;

           pairwiseBW.write("\t");
           for (int taxa = 0; taxa < theTaxa.length; taxa++) {

                pairwiseBW.write(theTaxa[taxa]+"\t");
            }
           pairwiseBW.write("\n");

           for (int j = 0; j < theTaxa.length; j++) {
               pairwiseBW.write(theTaxa[j]+"\t");

               for (int i = 0; i < theTaxa.length; i++) {
                   num= pairwiseMatrix[j][i];
                   pairwiseBW.write(Double.toString(num)+"\t");
               }
               pairwiseBW.write("\n");

           }

       pairwiseBW.close();
       }

       catch (IOException e) {
               System.out.println(e);
           }

       try {

            alignedPairwiseFile.createNewFile();
            alignedPairwiseFile.setWritable(true);
            BufferedWriter alignedPairwiseBW = new BufferedWriter(new FileWriter(alignedPairwiseFile));
            double numAligned = 0.0;

            alignedPairwiseBW.write("\t");
            for (int taxa = 0; taxa < theTaxa.length; taxa++) {

                alignedPairwiseBW.write(theTaxa[taxa]+"\t");
            }

            alignedPairwiseBW.write("\n");

            for (int j = 0; j < theTaxa.length; j++) {
                alignedPairwiseBW.write(theTaxa[j]+"\t");

               for (int i = 0; i < theTaxa.length; i++) {
                   numAligned= alignedPairwiseMatrix[j][i];
                   alignedPairwiseBW.write(Double.toString(numAligned)+"\t");
               }
               alignedPairwiseBW.write("\n");

           }

       alignedPairwiseBW.close();
       }

       catch (IOException e) {
               System.out.println(e);
           }
   }

/**Compares the sequence of tags in an input file to tags in master file (such as a build file) and outputs tags in the build file with one and two mismatches
 * in the 64bp sequence. Does not provide filters. Filter both input and master tag files before using. Because the master file is typically very large this
 * program is slow. For speed, use input tag file with only tags of interest (ie, those output from NovelTagsInTBTOutputToTBT). Input files can be either TBT
 * or TagCount. Output to TagCount.**/
   public static void NearMismatchInTagCount() {
       String inTBTFileName= "/Users/kelly/Documents/GBS/SW_landraces/output/"+compareFileID;
       String compareTagCountFileName= "/Users/kelly/Documents/GBS/landraceBuild/"+masterFileID;
//        String sameTagCountFileName= "/Users/kelly/Documents/GBS/landraceBuild/"+compareFileID+masterFileID+"TagCountIdentical";
       String oneTagCountFileName= "/Users/kelly/Documents/GBS/landraceBuild/"+compareFileID+masterFileID+"TagCountOneMismatch";
       String twoTagCountFileName= "/Users/kelly/Documents/GBS/landraceBuild/"+compareFileID+masterFileID+"TagCountTwoMismatch";
       File outFile= new File("/Users/kelly/Documents/GBS/landraceBuild/"+compareFileID+masterFileID+"TagCountOutput");
       TagCounts inTags= new TagCounts(inTBTFileName, FilePacking.Bit);
       TagCounts compareTags= new TagCounts(compareTagCountFileName, FilePacking.Bit);
//        TagCountMutable perfectMatch= new TagCountMutable(inTags, inTags.getTotalCount());
       TagCountMutable oneMismatch= new TagCountMutable(compareTags.getTag(0).length, compareTags.getTagCount());
       TagCountMutable twoMismatch= new TagCountMutable(compareTags.getTag(0).length, compareTags.getTagCount());
       String currTag;
       String compareTag;
       long identical= 0;
       int mismatch= 0;
       long oneNearMatch= 0;
       long twoNearMatch= 0;
       int tagsEval= 0;
       int read= 0;
       int[][] matching= new int[inTags.getTagCount()][3];//first holds perfect matches, second holds 1 mismatch, third holds 2 mismatches

       for (int inputTag= 0; inputTag < inTags.getTagCount(); inputTag++) {
           currTag= BaseEncoder.getSequenceFromLong(inTags.getTag(inputTag));
           tagsEval++;
           for (int masterTag= 0; masterTag < compareTags.getTagCount(); masterTag++) {
               compareTag= BaseEncoder.getSequenceFromLong(compareTags.getTag(masterTag));
               while (read < currTag.length() && mismatch<4) {
                   if (currTag.charAt(read)!=compareTag.charAt(read)) mismatch++;
                   read++;
               }
               if(mismatch == 0) {
   //                perfectMatch.addReadCount(currTag, currTag.length, inTags.getReadCount(inputTag));
                   matching[inputTag][0]++;
                   identical++;
               }
               else if(mismatch == 1) {
                   oneMismatch.addReadCount(compareTags.getTag(inputTag), compareTags.getTagSizeInLong(), compareTags.getReadCount(inputTag));
                   matching[inputTag][1]++;
                   oneNearMatch++;
               }
               else if(mismatch == 2) {
                   twoMismatch.addReadCount(compareTags.getTag(inputTag), compareTags.getTagSizeInLong(), compareTags.getReadCount(inputTag));
                   matching[inputTag][2]++;
                   twoNearMatch++;
               }
               mismatch= 0;
               read= 0;
           }
           if (tagsEval%10==0) System.out.println("TagsEval: "+tagsEval+"\tPerfectMatches: "+identical+"\tOneMismatch: "+oneNearMatch+"\tTwoMismatch: "+twoNearMatch);
       }
       System.out.println("TagsEvaluated: "+tagsEval+"\tPerfectMatches: "+identical+"\tOneMismatch: "+oneNearMatch+"\tTwoMismatch: "+twoNearMatch);

//        perfectMatch.writeTagCountFile(sameTagCountFileName, FilePacking.Bit, 1);
       oneMismatch.writeTagCountFile(oneTagCountFileName, FilePacking.Bit, 1);
       twoMismatch.writeTagCountFile(twoTagCountFileName, FilePacking.Bit, 1);

       try {
           outFile.createNewFile();
           outFile.setWritable(true);
           BufferedWriter BW = new BufferedWriter(new FileWriter(outFile));

           BW.write("TagSequence\tPerfectMatches\tOneMismatch\tTwoMismatches\n");
           for (int h=0; h < inTags.getTagCount(); h++) {
               BW.write(BaseEncoder.getSequenceFromLong(inTags.getTag(h))+"\t"+matching[h][0]+"\t"+matching[h][1]+"\t"+matching[h][2]+"\n");
           }
           BW.close();
       }

       catch (IOException e) {
           System.out.println(e);
       }

   }

   public static void main (String args[]) {
//        NovelTagsInTBT();
//        NovelTagsInTBTOutputToTBT(10,1);

//       int[] tags= new int[]{2,3,4,5};
//       int[] taxa= new int[]{1};
//       NovelTagsInTBTOutputToTBTReadsPerTag(tags,taxa);
       
       CDFOrderedTagCounts();

//        PairwiseTagComparisonsInTBT();
//        PairwiseTagComparisonsWithSharedReadCorrection();
//        SymmetricalPairwiseTagComparisonsWithSharedReadCorrection();
//        filterMergedTBT(3);
//        compareFileID= "AlignedNovelTagsFromSW_landraces_min1.tbt.byte_minReads5_minTaxa2";
//        masterFileID= "build20120110MergedTags.cnt";
//        NearMismatchInTagCount();
//        compareFileID= "AlignedNovelTagsFromSW_landraces_min1.tbt.byte_minReads5_minTaxa3";
//        masterFileID= "build20120110MergedTags.cnt";
//        NearMismatchInTagCount();
//        compareFileID= "AlignedNovelTagsFromSW_landraces_min1.tbt.byte_minReads10_minTaxa1";
//        masterFileID= "build20120110MergedTags.cnt";
//        NearMismatchInTagCount();
//        compareFileID= "AlignedNovelTagsFromSW_landraces_min1.tbt.byte_minReads5_minTaxa1";
//        masterFileID= "build20120110MergedTags.cnt";
//        NearMismatchInTagCount();
//        compareFileID= "UnalignedNovelTagsFromSW_landraces_min1.tbt.byte_minReads5_minTaxa2";
//        masterFileID= "build20120110MergedTags.cnt";
//        NearMismatchInTagCount();
//        compareFileID= "UnalignedNovelTagsFromSW_landraces_min1.tbt.byte_minReads5_minTaxa3";
//        masterFileID= "build20120110MergedTags.cnt";
//        NearMismatchInTagCount();
//        compareFileID= "UnalignedNovelTagsFromSW_landraces_min1.tbt.byte_minReads10_minTaxa1";
//        masterFileID= "build20120110MergedTags.cnt";
//        compareFileID= "UnalignedNovelTagsFromSW_landraces_min1.tbt.byte_minReads5_minTaxa1";
//        masterFileID= "build20120110MergedTags.cnt";
//        NearMismatchInTagCount();
   }

}