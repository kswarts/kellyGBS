/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package Kelly;

import java.io.BufferedReader;
import java.util.ArrayList;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.*;
import net.maizegenetics.util.Utils;

/**
 *
 * @author kelly
 */
public class BetterTagMethods {
    private static int numPop= 0;
    private static int[] popGroup;
    private static String dir;
    
    private static void filterUnmatchedTags(String tbtFileName, String topmFileName, String tagCountFileName, String pop, int minTotalRead, int minPopSpecific) {
       TagCounts theTags = new TagCounts(tagCountFileName, TagsByTaxa.FilePacking.Bit);//reads in the master tag file specified
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(topmFileName, true);
       TagsByTaxa theTBT= new TagsByTaxaByte(tbtFileName, TagsByTaxa.FilePacking.Byte);
       String[] taxaNames= theTBT.getTaxaNames();
       readPopFile(pop, taxaNames);
       ArrayList<Integer> maps= new ArrayList<Integer>();
       ArrayList<Integer> noMap= new ArrayList<Integer>();
       int badTags= 0;
        for (int tag = 0; tag < theTBT.getTagCount(); tag++) {
            int[] readsForPop= new int[numPop];
            boolean anyGoodPops= false;
            if (theTags.getReadCount(tag)<minTotalRead) {//total read count for tag from new fastq does not meet threshold
                //check to see if it meets the pop min
                for (int taxon = 0; taxon < popGroup.length; taxon++) {
                    readsForPop[popGroup[taxon]]+=theTBT.getReadCountForTagTaxon(tag, taxon);
                }
                for (int reads = 0; reads < readsForPop.length; reads++) {
                    if (readsForPop[reads]>minPopSpecific) anyGoodPops= true;
                }
                if (anyGoodPops==false) {
                    badTags++;
                    continue;
                }
            }
            if (theTOPM.getEndPosition(tag)<0) {//if it maps to B73
                maps.add(tag);
            }
            else {//if it doens't map to B73
                noMap.add(tag);
            }
        }
        maps.trimToSize();
        noMap.trimToSize();
        System.out.println("Tags that do not meet minimum threshold for total number of population specific reads: "+badTags+
                "/nNovel tags that do map: "+noMap.size()+"\nNovel tags that map: "+maps.size());
        printToTagCnt(theTags, maps, "NovelsTagsThatMap.cnt");
        printToTagCnt(theTags, noMap, "NovelsTagsThatDoNotMap.cnt");
    }
    
    private static void readPopFile(String pop, String[] taxaNames) {
        BufferedReader fileIn = null;
        try {
            fileIn = Utils.getBufferedReader(pop, 1000000);
            popGroup= new int[taxaNames.length];
            for (int i = 0; i < taxaNames.length; i++) {
                String currLine = fileIn.readLine();
                char one= currLine.charAt(currLine.lastIndexOf(taxaNames[i].charAt(taxaNames[i].length()-1))+2);
                if (currLine.indexOf(one)==currLine.length()-1) {
                    popGroup[i]= Character.getNumericValue(one);
                    if (popGroup[i]>numPop) numPop= popGroup[i];
                    continue;
                }
                char two= currLine.charAt(currLine.lastIndexOf(taxaNames[i].charAt(taxaNames[i].length()-1))+3);
                String dd= Character.toString(one)+Character.toString(two);
                popGroup[i]= Integer.valueOf(dd);
                if (popGroup[i]>numPop) numPop= popGroup[i];
            }
        }
        catch (Exception e) {
        }
    }
    
    private static void printToTagCnt(TagCounts theTags, ArrayList goodTags, String outFileName) {
        TagCounts newTags= new TagCounts(theTags.getTagSizeInLong(), goodTags.size());
        int newIndex= 0;
        for (int tag = 0; tag < theTags.getTagCount(); tag++) {
            if (goodTags.contains(tag)) {
                newTags.setTag(theTags.getTag(tag), theTags.getTagLength()[tag], theTags.getReadCount(tag), newIndex);
                newIndex++;
            }    
        }
        newTags.writeTagCountFile(dir+outFileName, TagsByTaxa.FilePacking.Byte, 1);
    }
    
    public static void main (String args[]) {
        dir= "/Users/kelly/Documents/GBS/Imputation/SmallFiles/";
        String tbtFile= dir+"";
        String tagCnt= dir+"";
        String topm= dir+"";
        String pop= dir+"";
        filterUnmatchedTags(tbtFile, topm, tagCnt, pop, 20, 5);
        
        
//        String[] taxaNames= new String[]{"asds","sfdsf","sgetgrg","ytjty","qetrwet","myumujf","gertert","ferfetgrs","turrtut","vsgrhrtu","eagertgse","yuiulou","aetesrts","liutlio","htdrthtrhd"
//,"aetrsetre","juyiyukjuyf","thtrdhyt","ersyrsyt","thrtrsyhyr"};
//        readPopFile(dir+"samplePop.txt",taxaNames);
   }
}
