package Kelly;

/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */



import net.maizegenetics.gbs.pipeline.*;
import java.io.File;
import net.maizegenetics.baseplugins.ExtractHapmapSubsetPlugin;
import net.maizegenetics.gbs.maps.TagsOnPhysicalMap;
import net.maizegenetics.gbs.tagdist.TagCounts;
import net.maizegenetics.gbs.tagdist.TagsByTaxa;
import net.maizegenetics.gbs.tagdist.TagsByTaxa.FilePacking;
import net.maizegenetics.gbs.tagdist.TagsByTaxaByte;
import net.maizegenetics.gbs.tagdist.TagsByTaxaUtils;
import net.maizegenetics.pal.alignment.ExportUtils;
import net.maizegenetics.pal.alignment.ImportUtils;
import net.maizegenetics.pal.alignment.MutableNucleotideAlignmentHDF5;
import net.maizegenetics.pipeline.TasselPipeline;
import net.maizegenetics.prefs.TasselPrefs;


/**
 *
 * @author jcg233
 */
public class KellyPipelinesGeneric {
//   public static String dir= "//home/local/MAIZE/kls283";
   public static String dir= "//home/local/MAIZE/kls283";
   public static void main(String[] args) {
       TasselPrefs.putAlignmentRetainRareAlleles(false);
//       runFindMergeHaplotypesPlugin();
//       runMinorWindowViterbiImputationPlugin();
       runExtractHapmapSubsetPlugin();
//        convertTextTagCountsToBinary();
//        convertBinaryTagCountsToText();
//        convertBinaryTBTToText();
//        convertTOPMToText();
//        tagCountsToFastQ();
//        runQseqToTagCountPlugin();
//        runFastqToTagCountPlugin();
//        runMergeMultipleTagCountPlugin();
//        runQseqToTBTPlugin();
//        runFastqToTBTPlugin();
//        runFastqToTBTPluginWithTaxaCount();
//        runMergeTagsByTaxaFilesPlugin();
//        printSumCountsInTBTByTaxa();
//        mergeTaxaInTBT();
//       convertSAMToTOPM();
//        filterTOPMWithPositions();
//        runTagsToSNPByAlignmentMTPlugin();
//        runMergeDuplicateSNPsPlugin();
//        runGBSHapMapFiltersPlugin();
//        runBiParentalErrorCorrection();
//        runMergeIdenticalTaxaPlugin();
//        runBinaryToTextPlugin();
//        filterMergeAllZea();
//        imputeAllZea();
//        runQseqToHapMapPlugin();
//        runRawReadsToHapMapPlugin();
//        analyzeRI_landraces_minCount3();
//        analyzeRI_landraces_minCount10();
//        outputTBTforNovelTagsInTBT();
//        analyzeD09FYACXX_5_PstIMaize();
//        analyzeB08AAABXX_1_PstI_IBM();
//        analyzeB08AAABXX_1_PstI_IBM_TBTByte();
//        analyzeIBM94ApeKI();
//        testTOPMExpandMaxVariants();
//        runTagsToSNPByAlignmentPlugin();
//        filterHapMapForTaxa();
//        getHapMapReport();
   }

   public static void runFindMergeHaplotypesPlugin() {
       String dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
//       String base= "AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1";
       String base= "AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_HomoSegForTaxaGreaterThan.1Het";
       String[] testArgs = new String[] {
            "-hmp",   dir+base+".hmp.txt.gz",
            "-o",     dir+base+"_HaplotypeMerge4k_s+.hmp.txt.gz",//Output file(s) must include 's+.' plus will be replace by segment (0..(~sites/hapSize)\n"
            "-oE",    dir+base+"_HaplotypeMerge4kError.txt",//Optional file to record site by sites errors as the haplotypes are developed\n"
            "-sC",    "10",//Start chromosome\n"
            "-eC",    "10",// End chromosome\n"
            "-mxDiv",  "0.01",//    Maximum divergence from founder haplotype\n"
            "-hapSize","4000",//    Preferred haplotype block size in sites\n"
            "-minPres", "500", //    Minimum number of present sites within input sequence to do the search\n"
            "-maxHap",  "3000",//    Maximum number of haplotypes per segment\n"
            "-maxOutMiss",  "0.4",//  Maximum frequency of missing data in the output haplotype"
       };
       String[] args = testArgs;
       FindMergeHaplotypesPlugin plugin = new FindMergeHaplotypesPlugin();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }
   
   public static void runMinorWindowViterbiImputationPlugin() {
//       String dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
       String dir= "//Users/kelly/Documents/GBS/Imputation/SmallFiles/";
       String base= "RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1";
       String donor= "AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_HomoSegForTaxaGreaterThan.1Het_HaplotypeMergeWithExtraLandrace_sX.hmp.txt.gz";
       
//       String[] minMnCnt= {"15","20","25","30"};
//       String[] mxInbErr= {".01",".02",".03",".04",".05"};
//       String[] mxHybErr= {".003",".004",".005",".008",".01"};
       String[] minMnCnt= {"20"};
       String[] mxInbErr= {".01"};
       String[] mxHybErr= {".003"};
       
       for (int i = 0; i < minMnCnt.length; i++) {
           for (int j = 0; j < mxInbErr.length; j++) {
               for (int k = 0; k < mxHybErr.length; k++) {
            String out= "_masked_HomoSegBlock200HetWithExtras8k.minMtCntBASELINE"+minMnCnt[i]+".mxInbErr"+mxInbErr[j]+".mxHybErr"+mxHybErr[k]+".c10";   
            String[] testArgs = new String[] {
                 "-hmp",   dir+base+"_masked.hmp.txt.gz", //Input HapMap file(s) 'c+' to denote variable chromosomes\n"
                 "-d",     dir+donor, //Donor haplotype files 'c+s+' to denote sections\n"
                 "-o",     dir+base+out, //Output HapMap file(s) 'c+' to denote variable chromosomes\n"
     //            "-d",     dir+"AllZeaGBS_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1_HaplotypeMerge_s+.hmp.txt.gz",
     //            "-o",     dir+base+"_defaultDonor8k.minMtCnt30.c+.imp.mhmp.h5", //Output HapMap file(s) 'c+' to denote variable chromosomes\n"
                 "-sC",    "10",  //Start chromosome
                 "-eC",    "10",  //End chromosome\n"
                 "-minMnCnt",    minMnCnt[i],  //Minimum number of minor alleles in the search window (or "+minMajorRatioToMinorCnt+"X major) default to 20
                 "-mxInbErr",    mxInbErr[j],    //Maximum inbred error rate\n"
                 "-mxHybErr",    mxHybErr[k],    //Maximum hybrid error rate\n"
     //            "-inbNNOff",    "8",    //Whether to use inbred NN (default:"+inbredNN+")\n"
     //            "-hybNNOff",    "8",    //Whether to use both the hybrid NN (default:"+hybridNN+")\n"
     //            "-mxDonH",    "8",  //Maximum number of donor hypotheses to be explored (default: "+maxDonorHypotheses+")\n"
                 "-mnTestSite",    "20",  //Minimum number of sites to test for NN IBD (default:"+minTestSites+")\n"
            };
            String[] args = testArgs;
            MinorWindowViterbiImputationPlugin plugin = new MinorWindowViterbiImputationPlugin();
            plugin.setParameters(args);
            plugin.performFunction(null);
//            MutableNucleotideAlignmentHDF5 outHDF5= MutableNucleotideAlignmentHDF5.getInstance(dir+base+out+".imp.hmp.h5");
//            ExportUtils.writeToHapmap(outHDF5, true, dir+base+out+".hmp.txt.gz", '\t', null);
            //for 1:300 internal mask
            ImputationAccuracy.RunTest(ImportUtils.readFromHapmap(dir+base+".hmp.txt.gz", null), 
                ImportUtils.readFromHapmap(dir+base+out+".hmp.txt", null),null, false, 300,.6,.01,.2,dir+base+out+"Accuracy.txt");
            //for mask against 55k
//            ImputationAccuracy.RunTest(ImportUtils.readFromHapmap(dir+"RIMMA_282_SNP55K_AGPv2_20100513__S45391.chr10_matchTo_RIMMA_282_v2.6_MERGEDUPSNPS_20130513_chr10subset__minCov0.1.hmp.txt.gz",null),
//                    ImportUtils.readFromHapmap(dir+base+out+".hmp.txt",null),
//                    ImportUtils.readFromHapmap(dir+base+".hmp.txt.gz",null),
//                    true, -1, .6, .01, .2, dir+base+out+".Accuracy55k.txt");
               }
           }
       }
   }
   
   public static void runExtractHapmapSubsetPlugin() {
        String dir = dir= "/home/local/MAIZE/kls283/GBS/Imputation/";
        String inHDF5 =  "AllZeaGBSv27";
        
        MutableNucleotideAlignmentHDF5 mnah5= MutableNucleotideAlignmentHDF5.getInstance(dir+inHDF5+".hmp.h5");
        ExportUtils.writeToHapmap(mnah5, true, dir+inHDF5+".hmp.txt.gz", '\t', null);
        
            String[] args = new String[]{
                "-h", dir+inHDF5+".hmp.txt.gz",
                "-o", dir+"AmesGBSv27.hmp.txt.gz",
                "-p", dir+"Ames(no EP or GEM).txt",
                "-a", "2" // at least 2 "alleles" (actually, genotypes) = polymorphic
            };
            ExtractHapmapSubsetPlugin plugin = new ExtractHapmapSubsetPlugin(null);
            plugin.setParameters(args);
            plugin.performFunction(null);
    }
   
   public static void convertTextTagCountsToBinary() {
       String textTagCountsFileS =   "C:/Users/jcg233/Documents/Bioinformatics/NextGen/HapMapV2/test_RandomPairedEndToTBT/FakeTagCounts.txt";
       String binaryTagCountsFileS = "C:/Users/jcg233/Documents/Bioinformatics/NextGen/HapMapV2/test_RandomPairedEndToTBT/FakeTagCounts.bin";
       TagCounts tc = new TagCounts(textTagCountsFileS, FilePacking.Text);
       tc.sort();
       tc.writeTagCountFile(binaryTagCountsFileS, FilePacking.Bit, 1);
   }

   public static void convertBinaryTagCountsToText() {
       String binaryTagCountsFileS = "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.cnt";
       String textTagCountsFileS =   "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.cnt.txt";
       TagCounts tc = new TagCounts(binaryTagCountsFileS, FilePacking.Bit);
       tc.writeTagCountFile(textTagCountsFileS, FilePacking.Text, 1);
   }

   public static void convertBinaryTBTToText() {
       String workDir = "/Users/kelly/Documents/GBS/";
       String binaryTBTFileS = workDir + "/SW_landraces/output/AlignedNovelTagsReadPerTagFromSW_landraces_min1_minReads30_minTaxa2.tbt.byte";
       String textTBTFileS =   workDir + "/SW_landraces/output/AlignedNovelTagsReadPerTagFromSW_landraces_min1_minReads30_minTaxa2TBT.txt";
       TagsByTaxa tbt = new TagsByTaxaByte(binaryTBTFileS, FilePacking.Byte);
       File textTBTFile = new File(textTBTFileS);
       tbt.writeDistFile(textTBTFile, FilePacking.Text, 0);
   }

   public static void convertTOPMToText() {
       String inTOPMFileName= "/Users/kelly/Documents/GBS/landraceCombo/topm/SW_RI_Span_landraces_min2.topm";
       String outTOPMFileName= "/Users/kelly/Documents/GBS/landraceCombo/topm/SW_RI_Span_landraces_min2.txt";
       File outTOPM= new File(outTOPMFileName);
       TagsOnPhysicalMap theTOPM= new TagsOnPhysicalMap(inTOPMFileName, true);
       theTOPM.writeTextFile(outTOPM);
   }

   public static void runFastqToTBTPluginWithTaxaCount() {

       //String testFastqDir = "H:/NAM_ApeKI_plates49_50/";
       String[] RI_landraces_minCount5_Args = new String[] {
           "-i", dir+"/GBS/B73/fastq",
           "-k", "/Users/kelly/Documents/GBS/RI_landraces/fastq/C08JYACXX_2_barcode_key.txt",
           "-e", "ApeKI",
           "-o", dir+"/GBS/B73/tbt",
           "-c", "1",
           "-y",
           "-t",dir+"/GBS/B73/mergedTagCounts/B73_min10.cnt",
       };

       String[] args = RI_landraces_minCount5_Args;
       FastqToTBTPluginWithTaxaCount plugin = new FastqToTBTPluginWithTaxaCount();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

  public static void runFastqToTBTPlugin() {

       //String testFastqDir = "H:/NAM_ApeKI_plates49_50/";
       String[] SW_landraces_minCount5_Args = new String[] {//run with minCount5 tagCount file and use short
           "-i", "/Users/kelly/Documents/GBS/fastq",
           "-k", "/Users/kelly/Documents/GBS/landraceCombo/fastq/Kelly_Landrace_key.txt",
           "-e", "ApeKI",
           "-o", "/Users/kelly/Documents/GBS/SW_landraces/tbt/w_Sorghum",
           "-c", "1",
           "-y",//outputs to tagsbyTaxaByte
           "-t","/Users/kelly/Documents/GBS/SW_landraces/mergedTagCounts/SW_landraces_w_SorghumMin1.cnt"
       };

       String[] RI_landraces_minCount1_Args = new String[] {//run with minCount1 tagCount file and use byte
           "-i", "/Users/kelly/Documents/GBS/RI_landraces/fastq",
           "-k", "/Users/kelly/Documents/GBS/RI_landraces/fastq/C08JYACXX_2_KLS_key.txt",
           "-e", "ApeKI",
           "-o", "/Users/kelly/Documents/GBS/RI_landraces/tbt",
           "-c", "1",//min taxa count to output
           "-y",//outputs to tagsbyTaxaByte
           "-t","/Users/kelly/Documents/GBS/RI_landraces/tagCounts/C08JYACXX_minCount1.cnt"
       };

       String[] Combo_landraces_minCount1_Args = new String[] {//run with minCount2 mergedTagCount file and use byte
           "-i", "/Users/kelly/Documents/GBS/fastq",
           "-k", "/Users/kelly/Documents/GBS/landraceCombo/fastq/Kelly_Landrace_key.txt",
           "-e", "ApeKI",
           "-o", "/Users/kelly/Documents/GBS/landraceCombo/tbt",
           "-c", "1",//min taxa count to output
           "-y",//outputs to tagsbyTaxaByte
           "-t","/Users/kelly/Documents/GBS/landraceCombo/mergedTagCounts/SW_RI_Span_landraces_min2.cnt"
       };
       
       String[] B73_Args = new String[] {
           "-i", dir+"/GBS/B73/fastq",
           "-k", dir+"/GBS/B73/fastq/B73__key_120611.txt",
           "-e", "ApeKI",
           "-o", dir+"/GBS/B73/tbt",
           "-c", "1",
           "-y",
           "-t",dir+"/GBS/B73/mergedTagCounts/B73_min10.cnt",
       };

       String[] args = B73_Args;
       FastqToTBTPlugin plugin = new FastqToTBTPlugin();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

   public static void runMergeTagsByTaxaFilesPlugin() {
       String[] SW_landraces = new String[] {
           "-i", "/Users/kelly/Documents/GBS/SW_landraces/tbt/w_Sorghum",
           "-o", "/Users/kelly/Documents/GBS/SW_landraces/mergedTBT/SW_landraces_wSorghum_min1.tbt.byte",
           "-x", //merges taxa with identical names
       };

       String[] Combo_landraces = new String[] {
           "-i", "/Users/kelly/Documents/GBS/landraceCombo/tbt",
           "-o", "/Users/kelly/Documents/GBS/landraceCombo/mergedTBT/SW_RI_Span_landrace_min1.tbt.byte",
           "-x", //merges taxa with identical names
       };

       String[] args = Combo_landraces;
       MergeTagsByTaxaFilesPlugin plugin = new MergeTagsByTaxaFilesPlugin();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

   public static void mergeTaxaInTBT() {
       String workDir = "/usr/local/maizediv/illumina/NAM_Ames_282/HapMap2MergedTBT/";
       String inputTBTFileS =            workDir + "HapMap2mergedTBT20110809.tbt.bin";
       String outputMergedTaxaTBTFileS = workDir + "HapMap2mergedTBT20110809.mergedTaxa.tbt.bin";
       TagsByTaxaUtils.mergeTaxaByName(inputTBTFileS, outputMergedTaxaTBTFileS, FilePacking.Bit, true);
       TagsByTaxaUtils.streamBinaryToText(outputMergedTaxaTBTFileS, 10000);
       TagsByTaxaUtils.printSumCounts(outputMergedTaxaTBTFileS, FilePacking.Bit, true);
   }

   public static void printSumCountsInTBTByTaxa() {
       String workDir = "/home/glaubitz/data/nextgen/allZea/tbt/";
       String inputTBTFileS = workDir + "allZea20110811.tbt.bin";
       TagsByTaxaUtils.printSumCounts(inputTBTFileS, FilePacking.Bit, true);
   }

   public static void convertSAMToTOPM() {
       
//        String SAMFile =      "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.sam";
//        String TOPMFile =     "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.topm.bin";
//        String TOPMTextFile = "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.topm.txt";

       String SAMFile =      dir+"/GBS/landraceBuild/AllZeaMasterTags_c10_20120607.sam";
       String TOPMFile =     dir+"/GBS/landraceBuild/AllZeaMasterTags_c10_20120607.topm";
//        String TOPMTextFile = "H:/70MU0AAXX/newestPipeline/NAM49_50_PstI_mergedTags_min50.topm.txt";

       TagsOnPhysicalMap topm = new TagsOnPhysicalMap();
       topm.readSAMFile(SAMFile, topm.getTagSizeInLong());
       topm.writeBinaryFile(new File(TOPMFile));
//        topm.writeTextFile(new File(TOPMTextFile));
   }
//
   public static void filterTOPMWithPositions() {
       String TOPMFile =     "N:/cassava/cassava.topm.bin";
       String TOPMFileWPos =     "N:/cassava/cassava_wPos.topm.bin";
       String TOPMTextFileWPos = "N:/cassava/cassava_wPos.topm.txt";

       TagsOnPhysicalMap topm = new TagsOnPhysicalMap(TOPMFile, true);
       topm.writeBinaryFile(new File(TOPMFileWPos), Integer.MAX_VALUE, true, true, (float) 1.0, true);
       topm = new TagsOnPhysicalMap(TOPMFileWPos, true);
       topm.writeTextFile(new File(TOPMTextFileWPos));
   }

   public static void runMergeDuplicateSNPsPlugin() {

       String testWorkdir = "C:/Users/jcg233/Documents/Bioinformatics/NextGen/NAM_ApeKI_plates49_50/hapmapOutput/testMergeDuplicateSNPsPlugin/";
       String[] testArgs = new String[] {
           "-hmp", testWorkdir+"mergedTBT_testSNPMerge_inputTestEnd_c+.hmp.txt", // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
           "-o",   testWorkdir+"mergedTBT_testSNPMerge_outputTestEnd_delUnmergDups_c+.hmp.txt", // Output HapMap file
           "-misMat", "0.15", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
           "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//          "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
           "-s", "10",        // Start chromosome (default 1)
           "-e", "10"         // End chromosome (default 10)
       };

       String workdir = "/Users/kelly/Documents/GBS/SW_landraces";
       String[] allZeaArgs = new String[] {
           "-hmp", workdir+    "rawReads2hapmap/.c+.hmp.txt", // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
           "-o",   workdir+"rawReads2hapmap/mergedSNPs/allZea_112311_SNPmerge15_c+.hmp.txt", // Output HapMap file
           "-misMat", "0.15", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
           "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//          "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
           "-s", "1",        // Start chromosome (default 1)
           "-e", "10"         // End chromosome (default 10)
       };

       String[] IBM94PstIArgs = new String[] {
           "-hmp",     "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/filtLD/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_c+.hmp.txt",
           "-o",       "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/filtLD/mergedSNPs/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_mergSNP20_c+.hmp.txt",
           "-misMat",  "0.2", // Threshold genotypic mismatch rate above which the duplicate SNPs won't be merged (default: 0.05)
           "-callHets",       // When two genotypes disagree at a SNP, call it a heterozygote (default: off = set to missing)
//            "-kpUnmergDups",   // When two duplicate SNPs were not merged (different alleles or too many mismatches), keep them (default: off = delete them)
           "-s", "1",        // Start chromosome (default 1)
           "-e", "10"         // End chromosome (default 10)
       };

       String[] args = testArgs;
       MergeDuplicateSNPsPlugin plugin = new MergeDuplicateSNPsPlugin();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

   public static void runGBSHapMapFiltersPlugin() {

       String workdir = "/usr/local/maizediv/illumina/allZeaHapMap/";

       String[] argsNAM49_50ApeKI = new String[] {
           "-hmp",    "H:/NAM_ApeKI_plates49_50/newestPipeline/hapmapOutput/20110726/mergedTBT.c+.hmp.txt", // Input HapMap file; use a plus sign (+) as a wild card character to specify multiple chromosome numbers
           "-o",      "H:/NAM_ApeKI_plates49_50/newestPipeline/hapmapOutput/20110726/filtered/NAM49_50_ApeKI_filtered.c+.hmp.txt", // Output HapMap file
           "-mnTCov", "0.05", // Minimum taxa coverage
           "-mnSCov", "0.05", // Minimum presence
           "-mnF",    "0.8",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
           "-mnMAF",  "0.1",  // Minimum minor allele frequency (default 0.0)
           "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
           "-hLD",            // Filter for high LD
           "-sC",     "1",    // Start chromosome (default 1)
           "-eC",     "10"    // End chromosome (default 10)"
       };

       String[] allZeaArgs = new String[] {
           "-hmp", workdir + "unfiltered/allZea20110812_unfiltSNPs_f9maf005mac20cov10_chr+.hmp.txt",
           "-o",   workdir +   "filtered/allZea20110812_filtF9maf001mac20LCov10_chr+.hmp.txt",
//            "-mnTCov", "0.0001",
           "-mnF",    "0.9",
           "-mnMAF",  "0.001",
//            "-hLD",
           "-mnSCov", "0.10", // Minimum locus coverage (proportion of Taxa)
           "-sC",     "8",    // Start chromosome
           "-eC",     "8"    // End chromosome
       };

       String[] BrunetArgs = new String[] {
           "-hmp", "C:/Users/jcg233/Documents/Bioinformatics/NextGen/UserSupport/Brunet/BC364mergedSNPs.c+.hmp.txt",
           "-o",   "C:/Users/jcg233/Documents/Bioinformatics/NextGen/UserSupport/Brunet/filt/BC364mergedSNPs_filt.c+.hmp.txt",
           "-mnTCov", "0.1",
//            "-mnF",    "0.9",
           "-mnMAF",  "0.0",
//            "-hLD",
           "-mnSCov", "0.1", // Minimum locus coverage (proportion of Taxa)
           "-sC",     "100",    // Start chromosome
           "-eC",     "101"    // End chromosome
       };

       String[] IBM94PstIArgs = new String[] {
           "-hmp",    "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinLCov10/filt/IBM94PstI_mnF80_filt_mnT10_mnS10_c+.hmp.txt",
           "-o",      "/usr/local/maizediv/illumina/Zea/PstI/hapmap/quantMnF80MinSCov80/filtLD/IBM94PstI_mnF80_filt_mnT10_mnS80_LD_c+.hmp.txt",
           "-mnTCov", "0.10", // Minimum taxa coverage (proportion of non-missing sites at a taxon)
           "-mnSCov", "0.80", // Minimum presence (proportion of non-missing taxa at a site)
           "-mnF",    "0.80",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
           "-mnMAF",  "0.20",  // Minimum minor allele frequency (default 0.0)
//            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
           "-hLD",            // Filter for high LD
           "-sC",     "1",    // Start chromosome (default 1)
           "-eC",     "10"    // End chromosome (default 10)
       };

       String[] IBM94ApeKIArgs = new String[] {
           "-hmp",    "/usr/local/maizediv/illumina/Zea/IBM94ApeKI/hapmap/filt/IBM94ApeKI_mnF80_filt_mnT03_mnS05_c+.hmp.txt",
           "-o",      "/usr/local/maizediv/illumina/Zea/IBM94ApeKI/hapmap/filtLD/IBM94ApeKI_mnF80_filt_mnT03_mnS05_LD_c+.hmp.txt",
           "-mnTCov", "0.03", // Minimum taxa coverage (proportion of non-missing sites at a taxon)
           "-mnSCov", "0.05", // Minimum presence (proportion of non-missing taxa at a site)
           "-mnF",    "0.80",  // Minimum F (inbreeding coefficient) (default -2.0  - no filter)
           "-mnMAF",  "0.20",  // Minimum minor allele frequency (default 0.0)
//            "-mxMAF",  "0.9",  // Maximum minor allele frequency (default 1.0 - no filter)
           "-hLD",            // Filter for high LD
           "-sC",     "1",    // Start chromosome (default 1)
           "-eC",     "10"    // End chromosome (default 10)
       };

       String[] args = IBM94ApeKIArgs;
       GBSHapMapFiltersPlugin plugin = new GBSHapMapFiltersPlugin();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

    public static void runMergeIdenticalTaxaPlugin() {
       String basedir = "/usr/local/maizediv/illumina/Zea/build20120110/bpec/TeoDNA_P1/";
       String[] TeoDNA_P1Args = new String[] {
           "-hmp", basedir+           "TeoDNA_P1_20120110_scv10mF8maf002_mgs_E1pLD5kpUn_chr+.hmp.txt",
           "-o",   basedir+"mergedTaxa/TeoDNA_P1_20120110_scv10mF8maf002_mgs_E1pLD5kpUn_mgt8_chr+.hmp.txt",
//            "-xHets",
           "-hetFreq", "0.8",
           "-sC", "1",  // Start chromosome
           "-eC", "10" // End chromosome
       };

       String[] args = TeoDNA_P1Args;
       MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
       mitp.setParameters(args);
       mitp.performFunction(null);
   }

   public static void runBinaryToTextPlugin() {
       String NAM_ApeKI_plates49_50Dir = "H:/NAM_ApeKI_plates49_50/newPipeline/";
       String[] TagCountsArgs = new String[] {
           "-i", NAM_ApeKI_plates49_50Dir+"mergedTagCounts/mergedApeKINAM49_50_min10.cnt",
           "-o", NAM_ApeKI_plates49_50Dir+"mergedTagCounts/mergedApeKINAM49_50_min10_testB2Tplugin.cnt.txt",
           "-t", "TagCounts",  // TOPM, TagCounts, TBTBit
       };

       String[] PstIMin5Args = new String[] {
           "-i", "/usr/local/maizediv/illumina/Zea/PstI/topm/B08AAABXX_1_IBM94PstI_min5.topm.bin",
           "-o", "/usr/local/maizediv/illumina/Zea/PstI/topm/B08AAABXX_1_IBM94PstI_min5.topm.txt",
           "-t", "TOPM",  // TOPM, TagCounts, TBTBit
       };

       String[] TestArgs = new String[] {
           "-i", "H:/NAM_ApeKI_plates49_50/testFastq/mergedTBT/testMerged2.tbt.bin",
           "-o", "H:/NAM_ApeKI_plates49_50/testFastq/mergedTBT/testMerged2.tbt.txt",
           "-t", "TBTBit",  // TOPM, TagCounts, TBTBit
       };

       String[] ProdTOPM20120110Args = new String[] {
           "-i", "/usr/local/maizediv/illumina/Zea/build20120110/topm/zea20120110c1-4.prod1-4.topm",
           "-o", "/usr/local/maizediv/illumina/Zea/build20120110/topm/zea20120110c1-4.prod1-4.topm.txt",
           "-t", "TOPM",  // TOPM, TagCounts, TBTBit
       };

       String[] args = ProdTOPM20120110Args;
       BinaryToTextPlugin plugin = new BinaryToTextPlugin(null);
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

   public static void filterMergeAllZea() {
       String workDir = "/usr/local/maizediv/illumina/Zea/hapmap/build111123/";
       String filePrefix = "allZea_111123_SNPmerge15";
       int sC=1;
       int eC=9;

       String[] args = new String[] {
           "-hmp", workDir+"mergedSNPs/"+filePrefix+"_c+.hmp.txt",
           "-o",       workDir+"filt/"+filePrefix+"_cov10_fT1_c+.hmp.txt",
           "-mnTCov", "0.01",
           "-mnF",    "0.9",
           "-mnMAF",  "0.002",
//            "-hLD",
           "-mnSCov", "0.10", // Minimum locus coverage (proportion of Taxa)
           "-sC", ""+sC,  // Start chromosome
           "-eC", ""+eC // End chromosome
       };
       GBSHapMapFiltersPlugin testClass = new GBSHapMapFiltersPlugin();
       testClass.setParameters(args);
       testClass.performFunction(null);

       args = new String[] {
           "-hmp", workDir+"filt/"+filePrefix+"_cov10_fT1_c+.hmp.txt",
           "-o",     workDir+"bpec/"+filePrefix+"_cov10_fT1E1pLD_c+.hmp.txt",
           "-oB", "usr/local/maizediv/illumina/allZeaHapMap/build111003/bpec/errorBin.txt",
           "-oE", "usr/local/maizediv/illumina/allZeaHapMap/build111003/bpec/errorBySNP.txt",
           "-popM", "Z[0-9]{3}",
           "-sC",    ""+sC,
           "-eC",    ""+eC,
           "-mxE",   "0.01",
           "-mnD",   "2.0",
           "-mnPLD", "0.5"
       };
       BiParentalErrorCorrectionPlugin bpec = new BiParentalErrorCorrectionPlugin();
       bpec.setParameters(args);
       bpec.performFunction(null);

       args = new String[] {
           "-hmp",      workDir+"bpec/"+filePrefix+"_cov10_fT1E1pLD_c+.hmp.txt",
           "-o",  workDir+"mergedTaxa/"+filePrefix+"_cov10_fT1E1pLD_mergedTaxa_c+.hmp.txt",
           "-xHets",
           "-hetFreq", "0.76",
           "-sC", ""+sC,  // Start chromosome
           "-eC", ""+eC // End chromosome
       };
       MergeIdenticalTaxaPlugin mitp = new MergeIdenticalTaxaPlugin();
       mitp.setParameters(args);
       mitp.performFunction(null);
   }

    public static void imputeAllZea() {
       int sC=2;
       int eC=2;
       String stemIn =  "/usr/local/maizediv/illumina/Zea/hapmap/build111123/mergedTaxa/allZea_112311_SNPmerge15";
       String stemOut = "/usr/local/maizediv/illumina/Zea/hapmap/build111123/imputed/allZea_112311_SNPmerge15";
       String[] args = new String[] {
           "-hmp", stemIn+"_cov10_fT1E1pLD_mergedTaxa_c+.hmp.txt",
           "-o", stemOut+"_cov10_fT1E1pLD_mergedTaxa_imputed_c+.hmp.txt",
           "-sC", ""+sC,  // Start chromosome
           "-eC", ""+eC // End chromosome
       };
       FastImputationBitFixedWindow.main(args);
   }

   public static void runRawReadsToHapMapPlugin() {

       String workDir = "/Users/kelly/Documents/GBS";
       String[] SW_landraces = new String[] { //run 20120319
           "-i", workDir + "/fastq/",   // contains a single HiSeq lane from DTMA (B08G4ABXX lane 1)
           "-k", workDir + "/SW_landraces/fastq/KLS15to16wSBkey.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
           "-o", workDir + "/SW_landraces/production",
           "-m", workDir + "/topm_filtered_042012.topm", // master TOPM file with variants recorded from discovery phase
//            "-d", "0",  // Maximum divergence between new read and previously mapped read (default = 0)
       };

       String[] RI_landraces = new String[] { //run 20120319
           "-i", workDir + "/fastq/",   // contains a single HiSeq lane from DTMA (B08G4ABXX lane 1)
           "-k", workDir + "/RI_landraces/fastq/C08JYACXX_2_KLS_key.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
           "-o", workDir + "/RI_landraces/",
           "-m", workDir + "/topm_filtered_042012.topm", // master TOPM file with variants recorded from discovery phase
//            "-d", "0",  // Maximum divergence between new read and previously mapped read (default = 0)
       };

       String[] NIL28Args = new String[] {
           "-i", "/usr/local/maizediv/illumina/Zea/NIL28/fastq",   // contains symbolic links to the 4 lanes of data
           "-k", "/usr/local/maizediv/illumina/Zea/NIL28/NIL28BarcodeKey.txt",
           "-e", "ApeKI",
           "-o", "/usr/local/maizediv/illumina/Zea/NIL28/ProdHapmap/oldTOPM",
           "-m", "/usr/local/maizediv/illumina/Zea/DTMA/topm/mergedNAM282Ames_variants2_20111018.topm.bin", // old TOPM with variants recorded from discovery phase (NOT COMPLEMENTED)
//            "-m", "/usr/local/maizediv/illumina/Zea/build20120110/topm/zea20120110c1-4.prod1-4.topm", // master TOPM file with variants recorded from discovery phase
//            "-d", "0",  // Maximum divergence between new read and previously mapped read (default = 0)
       };

       String[] args = RI_landraces;
       RawReadsToHapMapPlugin plugin = new RawReadsToHapMapPlugin();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

   public static void analyzeRI_landraces_minCount3() {
       String[] args;
       String baseDir = "/Users/kelly/Documents/GBS/RI_landraces/";

       String[] FastqToTagCountArgs = new String[] {
           "-i", baseDir+"fastq",
           "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
//            "-c", "3", // Minimum tag count (default is 1).
           "-o", baseDir+"20120310/tagCounts",
       };
       args = FastqToTagCountArgs;
       FastqToTagCountPlugin pluginFastqToTagCount = new FastqToTagCountPlugin();
       pluginFastqToTagCount.setParameters(args);
       pluginFastqToTagCount.performFunction(null);

       String[] FastqToTBTArgs = new String[] {
           "-i", baseDir+"fastq",
           "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
           "-o", baseDir+"20120310/tbt",
           "-c", "3", // Minimum tag count (default is 1).
//            "-sh", "",//outputs to a short.  doesn't seem to work
           "-y", "", //outputs in short format
           "-t", baseDir+"20120310/tbt/C08JYACXX_2_min3.cnt", // master Tags file (only one lane)
       };
       args = FastqToTBTArgs;
       FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
       pluginFastqToTBT.setParameters(args);
       pluginFastqToTBT.performFunction(null);

       String[] MergeMultipleTagCountArgs = new String[] {
           "-i", baseDir+"20120310/tagCounts", // Input directory containing .cnt files
           "-o", baseDir+"20120310/mergedTagCounts/C08JYACXX_2_minCount3.tbt.shrt", //  Output file name
           "-c", "3" // Minimum count of reads to be output (default 1)
//            "-m" // Specifies that reads should be output in fastq.
       };
       args = MergeMultipleTagCountArgs;
       MergeMultipleTagCountPlugin pluginMergeMultipleTagCount = new MergeMultipleTagCountPlugin();
       pluginMergeMultipleTagCount.setParameters(args);
       pluginMergeMultipleTagCount.performFunction(null);

       // next is BWA...
   }

       public static void analyzeRI_landraces_minCount10() {
       String[] args;
       String baseDir = "/Users/kelly/Documents/GBS/RI_landraces/";

//        String[] FastqToTagCountArgs = new String[] {
//            "-i", baseDir+"fastq",
//            "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
//            "-e", "ApeKI", // Enzyme used to create the GBS library
////            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
////            "-c", "3", // Minimum tag count (default is 1).
//            "-o", baseDir+"20120310/tagCounts",
//        };
//        args = FastqToTagCountArgs;
//        FastqToTagCountPlugin pluginFastqToTagCount = new FastqToTagCountPlugin();
//        pluginFastqToTagCount.setParameters(args);
//        pluginFastqToTagCount.performFunction(null); ***don't need this plugin because use the one from minCOunt3 (which has minCount of 1

       String[] MergeMultipleTagCountArgs = new String[] {
           "-i", baseDir+"20120310/tagCounts", // Input directory containing .cnt files
           "-o", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min10.cnt", //  Output file name
           "-c", "10" // Minimum count of reads to be output (default 1)
//            "-m" // Specifies that reads should be output in fastq.
       };
       args = MergeMultipleTagCountArgs;
       MergeMultipleTagCountPlugin pluginMergeMultipleTagCount = new MergeMultipleTagCountPlugin();
       pluginMergeMultipleTagCount.setParameters(args);
       pluginMergeMultipleTagCount.performFunction(null);

       String[] FastqToTBTArgs = new String[] {
           "-i", baseDir+"fastq",
           "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
           "-o", baseDir+"20120311/tbt",
           "-c", "1", // Minimum tag count (default is 1).
//            "-sh", "",//outputs to a short.  doesn't seem to work
           "-y", "", //outputs in byte format
           "-t", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min10.cnt", // master Tags file (only one lane)
       };
       args = FastqToTBTArgs;
       FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
       pluginFastqToTBT.setParameters(args);
       pluginFastqToTBT.performFunction(null);

       // next is BWA...
   }

       public static void outputTBTforNovelTagsInTBT() {
       String[] args;
       String baseDir = "/Users/kelly/Documents/GBS/RI_landraces/";

//        String[] MergeMultipleTagCountMin1 = new String[] {
//            "-i", baseDir+"20120310/tagCounts", // Input directory containing .cnt files
//            "-o", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min1.cnt", //  Output file name
//            "-c", "1" // Minimum count of reads to be output (default 1)
////            "-m" // Specifies that reads should be output in fastq.
//        };
//        args = MergeMultipleTagCountMin1;
//        MergeMultipleTagCountPlugin pluginMergeMultipleTagCountMin1 = new MergeMultipleTagCountPlugin();
//        pluginMergeMultipleTagCountMin1.setParameters(args);
//        pluginMergeMultipleTagCountMin1.performFunction(null);

       String[] MergeMultipleTagCountMin2 = new String[] {
           "-i", baseDir+"20120310/tagCounts", // Input directory containing .cnt files
           "-o", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min2.cnt", //  Output file name
           "-c", "2" // Minimum count of reads to be output (default 1)
//            "-m" // Specifies that reads should be output in fastq.
       };
       args = MergeMultipleTagCountMin2;
       MergeMultipleTagCountPlugin pluginMergeMultipleTagCountMin2 = new MergeMultipleTagCountPlugin();
       pluginMergeMultipleTagCountMin2.setParameters(args);
       pluginMergeMultipleTagCountMin2.performFunction(null);

//        String[] MergeMultipleTagCountMin3 = new String[] {
//            "-i", baseDir+"20120310/tagCounts", // Input directory containing .cnt files
//            "-o", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min3.cnt", //  Output file name
//            "-c", "3" // Minimum count of reads to be output (default 1)
////            "-m" // Specifies that reads should be output in fastq.
//        };
//
//        args = MergeMultipleTagCountMin3;
//        MergeMultipleTagCountPlugin pluginMergeMultipleTagCountMin3 = new MergeMultipleTagCountPlugin();
//        pluginMergeMultipleTagCountMin3.setParameters(args);
//        pluginMergeMultipleTagCountMin3.performFunction(null);

//        String[] MergeMultipleTagCountMin5 = new String[] {
//            "-i", baseDir+"20120310/tagCounts", // Input directory containing .cnt files
//            "-o", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min5.cnt", //  Output file name
//            "-c", "5" // Minimum count of reads to be output (default 1)
////            "-m" // Specifies that reads should be output in fastq.
//        };
//        args = MergeMultipleTagCountMin5;
//        MergeMultipleTagCountPlugin pluginMergeMultipleTagCountMin5 = new MergeMultipleTagCountPlugin();
//        pluginMergeMultipleTagCountMin5.setParameters(args);
//        pluginMergeMultipleTagCountMin5.performFunction(null);

//        String[] MergeMultipleTagCountMin10 = new String[] {
//            "-i", baseDir+"20120310/tagCounts", // Input directory containing .cnt files
//            "-o", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min10.cnt", //  Output file name
//            "-c", "10" // Minimum count of reads to be output (default 1)
////            "-m" // Specifies that reads should be output in fastq.
//        };
//        args = MergeMultipleTagCountMin10;
//        MergeMultipleTagCountPlugin pluginMergeMultipleTagCountMin10 = new MergeMultipleTagCountPlugin();
//        pluginMergeMultipleTagCountMin10.setParameters(args);
//        pluginMergeMultipleTagCountMin10.performFunction(null);

//        String[] FastqToTBTMin1 = new String[] {
//            "-i", baseDir+"fastq",
//            "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
//            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-o", baseDir+"20120311/tbt/min1",
//            "-c", "1", // Minimum tag count (default is 1).
////            "-sh", "",//outputs to a short.  doesn't seem to work
//            "-y", "", //outputs in byte format
//            "-t", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min1.cnt", // master Tags file (only one lane)
//        };
//        args = FastqToTBTMin1;
//        FastqToTBTPlugin pluginFastqToTBTMin1 = new FastqToTBTPlugin();
//        pluginFastqToTBTMin1.setParameters(args);
//        pluginFastqToTBTMin1.performFunction(null);

       String[] FastqToTBTMin2 = new String[] {
           "-i", baseDir+"fastq",
           "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
           "-o", baseDir+"20120311/tbt/min2",
           "-c", "1", // Minimum tag count (default is 1).
//            "-sh", "",//outputs to a short.  doesn't seem to work
           "-y", "", //outputs in byte format
           "-t", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min2.cnt", // master Tags file (only one lane)
       };
       args = FastqToTBTMin2;
       FastqToTBTPlugin pluginFastqToTBTMin2 = new FastqToTBTPlugin();
       pluginFastqToTBTMin2.setParameters(args);
       pluginFastqToTBTMin2.performFunction(null);

//        String[] FastqToTBTMin3 = new String[] {
//            "-i", baseDir+"fastq",
//            "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
//            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-o", baseDir+"20120311/tbt/min3",
//            "-c", "1", // Minimum tag count (default is 1).
////            "-sh", "",//outputs to a short.  doesn't seem to work
//            "-y", "", //outputs in byte format
//            "-t", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min3.cnt", // master Tags file (only one lane)
//        };
//        args = FastqToTBTMin3;
//        FastqToTBTPlugin pluginFastqToTBTMin3 = new FastqToTBTPlugin();
//        pluginFastqToTBTMin3.setParameters(args);
//        pluginFastqToTBTMin3.performFunction(null);

//        String[] FastqToTBTMin5 = new String[] {
//            "-i", baseDir+"fastq",
//            "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
//            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-o", baseDir+"20120311/tbt/min5",
//            "-c", "1", // Minimum tag count (default is 1).
////            "-sh", "",//outputs to a short.  doesn't seem to work
//            "-y", "", //outputs in byte format
//            "-t", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min5.cnt", // master Tags file (only one lane)
//        };
//        args = FastqToTBTMin5;
//        FastqToTBTPlugin pluginFastqToTBTMin5 = new FastqToTBTPlugin();
//        pluginFastqToTBTMin5.setParameters(args);
//        pluginFastqToTBTMin5.performFunction(null);

//        String[] FastqToTBTMin10 = new String[] {
//            "-i", baseDir+"fastq",
//            "-k", baseDir+"fastq/C08JYACXX_2_barcode_key.txt",
//            "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-o", baseDir+"20120311/tbt/min10",
//            "-c", "1", // Minimum tag count (default is 1).
////            "-sh", "",//outputs to a short.  doesn't seem to work
//            "-y", "", //outputs in byte format
//            "-t", baseDir+"20120311/mergedTagCounts/C08JYACXX_2_min10.cnt", // master Tags file (only one lane)
//        };
//        args = FastqToTBTMin10;
//        FastqToTBTPlugin pluginFastqToTBTMin10 = new FastqToTBTPlugin();
//        pluginFastqToTBTMin10.setParameters(args);
//        pluginFastqToTBTMin10.performFunction(null);

   }

   public static void analyzeB08AAABXX_1_PstI_IBM_TBTByte() {
       String baseDir = "/usr/local/maizediv/illumina/Zea/PstI/";

       String[] FastqToTBTArgs = new String[] {
           "-i", baseDir+"fastq",
           "-k", baseDir+"B08AAABXX_1_IBM94PstI_key.txt",
           "-e", "PstI", // Enzyme used to create the GBS library
           "-o", baseDir+"tbtByte",
           "-c", "1", // Minimum tag count (default is 1).
           "-y", // use TagsByTaxaByte
           "-t", baseDir+"tagCounts/B08AAABXX_1.cnt", // master Tags file (only one lane)
       };
       FastqToTBTPlugin pluginFastqToTBT = new FastqToTBTPlugin();
       pluginFastqToTBT.setParameters(FastqToTBTArgs);
       pluginFastqToTBT.performFunction(null);

       // convery the binary tbt.byte file to a text version
       String tbtByteFile = baseDir+"tbtByte/B08AAABXX_1.tbt.byte";
       TagsByTaxa tbtByteIBMPstI = new TagsByTaxaByte(tbtByteFile, FilePacking.Byte);
       String tbtByteTextFile = baseDir+"tbtByte/B08AAABXX_1.tbt.txt";
       tbtByteIBMPstI.writeDistFile(new File(tbtByteTextFile), FilePacking.Text, 1);

       String[] TagsToSNPByAlignmentArgs = new String[] {
           "-i",    baseDir+"tbtByte/B08AAABXX_1.tbt.byte",
           "-y", // use TagsByTaxaByte
           "-o",    baseDir+"hapmap/quantMnF80MinLCov10",
           "-m",    baseDir+"topm/B08AAABXX_1_IBM94PstI_min50.topm.bin",
//            "-mUpd", baseDir+"",
           "-mnF",   "0.8",  // allow some more hets in IBM (tried 0.8 initially, then 0.6, 0.7, and 0.8 again)
           "-mnMAF", "0.2",
           "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
           "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
           "-s", "1",  // Start chromosome
           "-e", "10"  // End chromosome
       };
       String[] args = TagsToSNPByAlignmentArgs;
       TagsToSNPByAlignmentPlugin pluginTagsToSNPByAlignment = new TagsToSNPByAlignmentPlugin();
       pluginTagsToSNPByAlignment.setParameters(args);
       pluginTagsToSNPByAlignment.performFunction(null);
   }

   public static void analyzeIBM94ApeKI() {
       String baseDir = "/usr/local/maizediv/illumina/Zea/IBM94ApeKI/";

       String[] QseqToTagCountArgs = new String[] {
           "-i", baseDir+"qseq",
           "-k", baseDir+"C05F2ACXX_5_key.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
//            "-s", "200000000", // Max good reads per lane. (Optional. Default is 200,000,000)
           "-c", "5", // Minimum tag count (default is 1).
           "-o", baseDir+"tagCounts",
       };
       QseqToTagCountPlugin pluginQseqToTagCount = new QseqToTagCountPlugin();
       pluginQseqToTagCount.setParameters(QseqToTagCountArgs);
       pluginQseqToTagCount.performFunction(null);

       String[] QseqToTBTArgs = new String[] {
           "-i", baseDir+"qseq",
           "-k", baseDir+"C05F2ACXX_5_key.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
           "-o", baseDir+"tbt",
           "-c", "1", // Minimum tag count (default is 1).
           "-y", // use TagsByTaxaByte
           "-t", baseDir+"tagCounts/C05F2ACXX_5.cnt", // master Tags file (only one lane)
       };
       QseqToTBTPlugin pluginQseqToTBT = new QseqToTBTPlugin();
       pluginQseqToTBT.setParameters(QseqToTBTArgs);
       pluginQseqToTBT.performFunction(null);

       String TagCountFileName = baseDir+"tagCounts/C05F2ACXX_5.cnt";
       String FastQFileName    = baseDir+"tagCounts/C05F2ACXX_5_IBM94ApeKI_min5.fastq";
       TagCounts tc = new TagCounts();
       tc.toFASTQ(TagCountFileName, FastQFileName);

       // Comment out the next two steps until BWA has been run
       // At that point, comment out the above steps

       String[] SAMConverterPluginArgs = new String[] {
           "-i", baseDir+"tagCounts/C05F2ACXX_5_IBM94ApeKI_min5.sam",
           "-o", baseDir+     "topm/C05F2ACXX_5_IBM94ApeKI_min5.topm.bin",
       };
       SAMConverterPlugin pluginSAMConverter = new SAMConverterPlugin();
       pluginSAMConverter.setParameters(SAMConverterPluginArgs);
       pluginSAMConverter.performFunction(null);

       String[] TagsToSNPByAlignmentArgs = new String[] {
           "-i",    baseDir+"tbt/C05F2ACXX_5.tbt.byte",
           "-y", // use TagsByTaxaByte
           "-o",    baseDir+"hapmap/unfilt",
           "-m",    baseDir+"topm/C05F2ACXX_5_IBM94ApeKI_min5.topm.bin",
//            "-mUpd", baseDir+"",
           "-mnF",   "0.8",  // allow some more hets in IBM
           "-mnMAF", "0.2",
           "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
           "-mnLCov","0.05", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
           "-s", "1",  // Start chromosome
           "-e", "10"  // End chromosome
       };
       TagsToSNPByAlignmentPlugin pluginTagsToSNPByAlignment = new TagsToSNPByAlignmentPlugin();
       pluginTagsToSNPByAlignment.setParameters(TagsToSNPByAlignmentArgs);
       pluginTagsToSNPByAlignment.performFunction(null);
   }

   public static void  testTOPMExpandMaxVariants() {
       String TOPMwVariantsFile = "/usr/local/maizediv/illumina/Zea/DTMA/topm/mergedNAM282Ames_variants2_20111018.topm.bin";
       TagsOnPhysicalMap topm = new TagsOnPhysicalMap(TOPMwVariantsFile, true);
       topm.printRows(1000, true, true);
       topm.expandMaxVariants(8);
       topm.printRows(1000, true, true);
   }

   public static void runTagsToSNPByAlignmentPlugin() {

       String workdir= "/Users/kelly/Documents/GBS/";
       String[] SWLandracesArgs = new String[] {

           "-i",    workdir+"SW_landraces/mergedTBT/SW_landraces_wSorghum_min1.tbt.byte",
           "-o",    workdir+"SW_landraces/hapmap/wSB/unfilt",
           "-m",    workdir+"topm_filtered_042012.topm",
           "-y",
//            "-mUpd", workdir+"topm_filtered_04012_wSW3.topm",
//            "-mnF",   "0.9",
//            "-mnMAF", "0.02",
//            "-mnMAC", "10",
//            "-mnLCov","0.02", // Minimum locus coverage (proportion of Taxa)
//            "-inclRare",
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
           "-s", "1",  // Start chromosome
           "-e", "10"  // End chromosome
       };

       String baseDirIBM94PstI = "/usr/local/maizediv/illumina/Zea/PstI/";
       String[] IBM94PstIArgs = new String[] {
           "-i",    baseDirIBM94PstI+"tbtByte/B08AAABXX_1.tbt.byte",
           "-y", // use TagsByTaxaByte
           "-o",    baseDirIBM94PstI+"hapmap/quantMnF80MinLCov10",
           "-m",    baseDirIBM94PstI+"topm/B08AAABXX_1_IBM94PstI_min50.topm.bin",
//            "-mUpd", baseDir+"",
           "-mnF",   "0.8",  // allow some more hets in IBM (tried 0.8 initially, then 0.6, 0.7, and 0.8 again)
           "-mnMAF", "0.2",
           "-mnMAC", "99999",  // this will never be satified: this way -mnMAF overrides it
           "-mnLCov","0.1", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
//            "-callBiSNPsWGap",  //call sites with a biallelic SNP plus a gap (e.g., A/C/-)
           "-s", "1",  // Start chromosome
           "-e", "10"  // End chromosome
       };

       String baseDirCassava = "N:/cassava/";
       String[] CassavaArgs = new String[] {
           "-i",    baseDirCassava+"D09WGACXX_7.tbt.byte",
           "-y", // use TagsByTaxaByte
           "-o",    baseDirCassava+"hapmap/unfilt",
           "-m",    baseDirCassava+"cassava.topm.bin",
//            "-mUpd", baseDir+"",
//            "-mnF",   "-0.5",  // cassava is outbred
           "-mnMAF", "0.05",
           "-mnMAC", "9999",  // this will never be satified: this way -mnMAF overrides it
           "-mnLCov","0.40", // Minimum locus coverage (proportion of Taxa)
//            "-inclGaps",  // Include sites where major or minor allele is a GAP
           "-s", "1",  // Start chromosome
           "-e", "13000"  // End chromosome
       };

       String[] args = SWLandracesArgs;
       TagsToSNPByAlignmentPlugin plugin = new TagsToSNPByAlignmentPlugin();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

   public static void runQuantPipeline() {  // this method is not finished
       String[] NAM49_50_ApeKIargs = new String[] {
           "-i", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/qseq",
           "-k", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/NAM49_50_ApeKI_key.txt",
           "-e", "ApeKI", // Enzyme used to create the GBS library
           "-o", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/tbt",
           "-c", "1", // Minimum tag count (default is 1).
           "-t", "/cbsufsrv4/data1/maizediv/illumina/NAM_ApeKI_plates49_50/mergedTagCounts/NAM49_50_ApeKI_mergedTags_min10.cnt", // master Tags file
       };
       String[] args = NAM49_50_ApeKIargs;
       QseqToTBTPlugin plugin = new QseqToTBTPlugin();
       plugin.setParameters(args);
       plugin.performFunction(null);
   }

    public static void filterHapMapForTaxa() {
       int sC=1;
       int eC=5;
       String taxaListFileName = "/usr/local/maizediv/illumina/Zea/build20120110/TeoDNA_P1_taxaFullNames.txt";
       String baseDir          = "/usr/local/maizediv/illumina/Zea/build20120110/imp/";
       String infileS =  baseDir+"Zea20120110_scv10mF8maf002_mgs_E1pLD5kpUn.imp95_1024.c+.hmp.txt";
       String outfileS = baseDir+"TeoDNA_P1/TeoDNA_P1_20120110_scv10mF8maf002_mgs_E1pLD5kpUn_imp95_1024_chr+.hmp.txt";
       String infile, outfile;
       TasselPipeline tp = null;
       String[] args;
       for (int chr=sC; chr<=eC; ++chr) {
           infile=infileS.replace("+", ""+chr);
           outfile=outfileS.replace("+", ""+chr);
           args = new String[] {
               "-fork1",
               "-h", infile,
               "-includeTaxaInFile", taxaListFileName,
               "-taxaJoinStrict", "true",
               "-export", outfile,
               "-exportType", "Hapmap",
               "-runfork1",
           };
           tp = new TasselPipeline(args, null);
       }
   }
}