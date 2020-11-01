package expressionmatrix;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Class which represents an experiment Summary table using PatMaN Alignments.
 */
public class SummaryTable {

    private final ArrayList<SampleSummary> summary = new ArrayList<>();
    private final File[] srnaFasta;
    private final File[] patmanAlignment;
    private final String[] sampleName;
    private final String[] condition;
    private final String[] replicate;
    private final String outputPrefix;
    private final HashMap<String, Integer[]> masterExpressionTable;
    
    /**
     * Constructs and initializes a new SummaryTable.
     * @param srnaFasta ordered list of sRNA FASTA files - ordered by sample table.
     * @param patmanAlign ordered list of patmanAlignment files - ordered by sample table.
     * @param sampleName ordered list of sample names - ordered by sample table.
     * @param condition ordered list of conditions - ordered by sample table.
     * @param replicate ordered list of replicate IDs - ordered by sample table.
     * @param outputPrefix The output path and prefix of the experiment summary table.
     */
    public SummaryTable(File[] srnaFasta, File[] patmanAlign, String[] sampleName, 
            String[] condition, String[] replicate, String outputPrefix) {
        this.srnaFasta = srnaFasta;
        this.patmanAlignment = patmanAlign; 
        this.sampleName = sampleName;
        this.condition = condition;
        this.replicate = replicate;
        this.outputPrefix = outputPrefix;
        this.masterExpressionTable = new HashMap<>();
    }
    
    /**
     * Generates and experiment summary and outputs it to disk.
     */
    public void GenerateExperimentSummary(){
        for(int sampleIndex = 0; sampleIndex < sampleName.length; sampleIndex++) {
            System.out.println("Constructing Sample Summary for: "+sampleName[sampleIndex]);
            
            SampleSummary s = new SampleSummary(srnaFasta[sampleIndex], patmanAlignment[sampleIndex], 
                    sampleName[sampleIndex], condition[sampleIndex], replicate[sampleIndex]);
            s.countTotalSrnaReads();
            s.calculateMappingRate(); 
            s.generateSampleExpressionTable();
            
            //Iterate through the expression table for the sample and add it to the master expression table.
            HashMap<String, Integer> sampleExpressionTable = s.getSampeReferenceExpressionTable();
            //Get the keys for the sample's expression table.
            Iterator<String> itr = sampleExpressionTable.keySet().iterator();
            //While we have keys to look at.
            while(itr.hasNext()){
                //Get the next srnaSequence_key in the list.
                String mir_id = itr.next();
                //Get the mapped abundance for the mir for this sample.
                Integer mappedAbundance = sampleExpressionTable.get(mir_id);
                //If the master table already has this id - from the same id in another sample.
                if(this.masterExpressionTable.containsKey(mir_id)){
                    //then put the mapped abundance for this mir_id into the table at index sampleIndex (sampleIndex.e. the column for this sample.)
                    this.masterExpressionTable.get(mir_id)[sampleIndex] = mappedAbundance;
                }else{
                    //Else, we have not seen this mir_id before, so we need to make a new structure for it (sampleIndex.e. a new array which will hold the mapped abundance for each sample).
                    Integer[] mappedAbundanceList = new Integer[this.sampleName.length];
                    //Initialize each of the values to 0 - do this because we do knot know at run time at what point we are adding the new mir_id and all must be set to 0 if not previously seen.
                    for(int initalizeIndex = 0; initalizeIndex < this.sampleName.length; initalizeIndex++){ 
                        mappedAbundanceList[initalizeIndex] = 0; 
                    }
                    //Place the mapped abundance into the correct column for the sample.
                    mappedAbundanceList[sampleIndex] = mappedAbundance;
                    //Now put the mapped abundance list into the master reference expression table.
                    this.masterExpressionTable.put(mir_id, mappedAbundanceList);
                }
            }
            //Add the summary we just made to the list of summarys.
            this.summary.add(s);
        }
        //All done building the sample summarys, so out put the experiment summary.
        this.outputExperimentSummary();
    }
    
    /**
     * Output the Experiment summary table to file. 
     */
    private void outputExperimentSummary() {
        //Open a writer stream.
        try (BufferedWriter w = new BufferedWriter(new FileWriter(new File(this.outputPrefix+".ExperimentSummary.txt")))) {
            w.write("EXPERIMENT SUMMARY:");
            w.newLine();
            //Write the sample name followed by the sample names in order.
            w.write("Sample name:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getSampleName()+"\t");
            }
            w.newLine();
            //Write the conition followed by the conditions in order.
            w.write("Condition:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getCondition()+"\t");
            }
            w.newLine();
            //Write the iDEP replicated id for each sample in order.
            w.write("Replicate_ID:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getReplicate()+"\t");
            }
            w.newLine();
            //Write the sRNA file name for each sample in order.
            w.write("SRNA_FASTA:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.SRNA_FASTA.getName()+"\t");
            }
            w.newLine();
            //Write the patman alignments file name for each sample in order.
            w.write("Alignments file name:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getPatmanAlign().getName()+"\t");
            }
            w.newLine();
            //Write the total srna count in each samples fasta file in sample order.
            w.write("Total read count in FASTA:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getTotalSrnasInSample()+"\t");
            }
            w.newLine();
            //write the number of srna reads which aligned in sample order.
            w.newLine();
            w.write("Number of reads aligned to the reference sequences:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getTotalSrnasWhichHaveMapping()+"\t");
            }
            w.newLine();
            //Give the break down of mapped once, twice (or more) or did not map in sample order.
            w.write("Of which:"); w.newLine();
            w.write("Mapped exactly once:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getTotalSrnasWhichHaveExactlyOneMapping()+"\t");
            }
            w.newLine();
            w.write("Mapped twice or more (multi-mappers):\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getTotalStransWichHaveTwoOrMoreMappings()+"\t");
            }
            w.newLine();
            w.newLine();
            w.write("Zero alignments found:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getTotalSrnasThatDidNotMap()+"\t");
            }
            w.newLine();
            //Give the break down of mapped once, twice (or more) or did not map in sample order as PERCENTAGES.
            w.newLine();
            w.write("Number of reads aligned to the reference sequences (%):\t");
            for (SampleSummary s : this.summary) {
               w.write(s.getTotalSrnasWhichHaveMapping_PERCENT()+"\t");
            }
            w.newLine();
            w.write("Of which:"); w.newLine();
            w.write("Mapped exactly once: (%):\t");
            for (SampleSummary s : this.summary) {
               w.write(s.getTotalSrnasWhichHaveExactlyOneMapping_PERCENT()+"\t");
            }
            w.newLine();
            w.write("Mapped twice or more (multi-mappers) (%):\t");
            for (SampleSummary s : this.summary) {
               w.write(s.getTotalStransWichHaveTwoOrMoreMappings_PERCENT()+"\t");
            }
            w.newLine();
            w.newLine();
            w.write("Zero alignments found (%):\t");
            for (SampleSummary s : this.summary) {
               w.write(s.getTotalSrnasWithoutMapping_PERCENT()+"\t");
            }
            w.newLine(); 
            
            //Give the number of reference sequences which had a mapping of one or more sRNAs.
            w.newLine(); 
            w.write("Number of reference sequences with read(s) aligned:\t");
            for (SampleSummary s : this.summary) {
               w.write(s.getExpressedReferenceCount()+"\t");
            }
            w.newLine(); 
            w.newLine();
            //Give the master reference expression table a title.
            w.write("REFERENCE SEQUENCE COUNTS:");
            w.newLine();
            //Print out the replicate id again - for the purpose of using with iDEP.
            w.write("Replicate_ID:\t");
            for (SampleSummary s : this.summary) {
                w.write(s.getReplicate()+"\t");
            }
            w.newLine();
            //Iterate over the master reference expression table and print it out.
            Iterator<String> itr = this.masterExpressionTable.keySet().iterator();
            while(itr.hasNext()){
                String mir = itr.next();
                Integer[] count = this.masterExpressionTable.get(mir);
                w.write(mir+"\t");
                for (Integer c : count) {
                    w.write(c + "\t");
                }
                w.newLine();
            }
        } catch (IOException ex) {
            System.err.println("Error: Cannot output the experiment summary to disk.");
            Logger.getLogger(SummaryTable.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
     * Class which represents a single sample within an experiment.
     */
    private class SampleSummary{
        /** Files and labels. **/
        private final File SRNA_FASTA;
        private final File patmanAlign;
        private final String sampleName;
        private final String condition;
        private final String replicate;
        /** The total number of sRNA reads in the sample.**/
        private Integer totalSrnasInSample;
        //COUNTS
        /** Total number of sRNA reads which have mapped. **/
        private Integer totalSrnasReadsWhichHaveMapping;
        /** Total number of sRNA reads which mapped exactly once. **/
        private Integer totalSrnaReadsWhichHaveExactlyOneMapping;
        /** Total number of sRNA reads which mapped twice or more. **/
        private Integer totalStranReadsWichHaveTwoOrMoreMappings;
        /** Total number of sRNA reads which did not map at all. **/
        private Integer totalSrnaReadsWithoutMapping;
        //PERCENTs
        /** Total number of sRNA reads which have mapped - as percent. **/ 
        private double totalSrnaReadsWhichHaveMapping_PERCENT;
        /** Total number of sRNA reads which mapped exactly once - as percent. **/
        private double totalSrnaReadsWhichHaveExactlyOneMapping_PERCENT;
        /** Total number of sRNA reads which mapped twice or more - as percent. **/
        private double totalStranReadsWichHaveTwoOrMoreMappings_PERCENT;
        /** Total number of sRNA reads which did not map at all - as percent. **/
        private double totalSrnaReadsWithoutMapping_PERCENT;
        /** The number of reference sequences which have one or more reads mapped to it. (Expressed reference sequences). **/
        private int expressedReference = 0;
        /** Total number of NR sequences which mapped once. **/
        private int mappedOnce = 0;
        /** Total number of NR sequences which were multi-mappers. **/
        private int multimapper = 0;   
        /** Reference expression table for this sample. **/
        HashMap<String, Integer> referenceExpressionMap = new HashMap<>();
        /**
         * ProcessedSrnaList = sRNAs that should not have any processing done as they have already been assigned
         * due to being a multi-mapper and the best alignments have already been assigned (i.e. those with 0 mismatches).
         */
        private final HashMap<String, String> processedSrnaList = new HashMap<>();
        
         
        /**
         * Constructs an instance of Sample Summary using supplied parameters.
         * @param SRNA_FASTA Small RNAs in NR Fasta format >sequence(abundance).
         * @param patmanAlign Patman alignment file which is not RC.
         * @param sampleName The name of the sample for which this sample summary represents.
         * @param condition The condition/treatment lable for this sample.
         * @param replicate The replicate ID for this sample (Important for iDEP/DEseq2).
         */
        public SampleSummary(File SRNA_FASTA, File patmanAlign, String sampleName, String condition, String replicate){
            this.SRNA_FASTA = SRNA_FASTA;
            this.patmanAlign = patmanAlign;
            this.sampleName = sampleName;
            this.condition = condition;
            this.replicate = replicate;
        }
         
        /**
         * Method counts the number of sRNAs reads within the NR FASTA file. 
         */
        private void countTotalSrnaReads() {
            //Set the total srnas in sample to zero.
            this.totalSrnasInSample = 0;
            //Open a reader stream to the fasta file for this sample.
            try (BufferedReader r = new BufferedReader(new FileReader(this.SRNA_FASTA))) {
                //Read the first alignment.
                String line = r.readLine();
                //While there are lines to read.
                while(line != null){
                    //Split the header alignment to get the abundance of the NR sRNA.
                    String[] toks = line.split("\\(");
                    String countS = toks[1].replace(")", "");
                    Integer abundance = Integer.parseInt(countS);
                    //Increment the total srna count by the abundance.
                    this.totalSrnasInSample += abundance;
                    //Read the following sequence alignment and discard it.
                    r.readLine();
                    //Read the next header alignment which contains the sequence and the abundance.
                    line = r.readLine();  
                }
            } catch (IOException ex) {
                System.out.println("Error: Could not count total srna reads within NR fasta file.");
                Logger.getLogger(SummaryTable.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
        
        /**
         * Method to build the expression matrix for the reference miRNAs from the PatMaN file provided to this class.
         * @return 
         */
        private void generateSampleExpressionTable() {
            //A place to store abundance counts for multi-mappers without double counting i.e. if there is an sRNA entry, then it is a multi-mapper.
            HashMap<String, Integer> multiMapperAbundance = new HashMap<>();
            //Map of <miRNA_ID:ReadCountsMappingToMiRNA>
            HashMap<String, Integer> referenceMap = new HashMap<>();
            //Open the patman file.
            try (BufferedReader r = new BufferedReader(new FileReader(this.patmanAlign))) {
                //Read the first alignment.
                String line = r.readLine();
                //While we have alignments to process.
                while(line != null){
                    //Split the current alignment into fields by tab.
                    String[] toks = line.split("\t");
                    //Get the mir_id ID.
                    String mir = toks[0];
                    //Get the sRNA and clean off the abundance - (note that patman file does not contain the > from the fasta).
                    String sRNA = toks[1].split("\\(")[0];
                    //Get the abundance for the sRNA mapping to the mir_id for this alignment.
                    Integer abundance = Integer.parseInt( toks[1].split("\\(")[1].replace(")", "") );
                    //Get the number of mismatches for this alignment.
                    Integer missmatches = Integer.parseInt(toks[5]);
                    //Query the mapping status for this sRNA - 0 is an exact mapper on a single miR, 1 is a multi-mapper to be checked for best .
                    int multimapperStatus = this.checkMultimapperStatus(sRNA);
                    //If the sRNA is NOT a multi-mapper
                    if(multimapperStatus == 0) {
                        if(referenceMap.containsKey(mir)){
                            Integer count = referenceMap.get(mir);
                            referenceMap.put(mir, count+abundance);
                        }else{
                            referenceMap.put(mir, abundance);
                            //First time we have seen a mapping to this mir_id, so increment the number of expressed mirs..
                            this.expressedReference++;
                        }
                        mappedOnce += abundance;
                    }
                    //It is a multi-mapper and we need to look to see if we can extract a best alignment or to allow all or to assign one (multi-map to single reference scenario).
                    else{
                        //Now to decide which alignment to keep.
                        //Check to see if the sRNA in the alignment has been previously identified as a multi mapper by checking in the processed sRNA list
                        //If already asgined do nothing as the best matches for the sRNA has already been processed.
                        if(!this.processedSrnaList.containsKey(sRNA)){
                            //Call the function findBestAlignments and add get the 'best matches' for this sRNA.
                            ArrayList<String> bestAlignments = this.findBestAlignments(line);
                            //If we have some best alignments, then add them to the reference map.
                            if(bestAlignments != null){
                                        
                                for(String align : bestAlignments){

                                    //Split the current best alignment into fields by tab.
                                    String[] toks_b = align.split("\t");
                                    //Get the current best mir_id.
                                    String mir_b = toks_b[0];
                                    //Get the sRNA for the current best alignment and clean off the abundance - (note that patman file does not contain the > from the fasta).
                                    String sRNA_b = toks_b[1].split("\\(")[0];
                                    //Get the abundance for the sRNA mapping to the mir_id for this current best alignment.
                                    Integer abundance_b = Integer.parseInt( toks_b[1].split("\\(")[1].replace(")", "") );
                                    //Get the number of mismatches for this current best alignment.
                                    Integer missmatches_b = Integer.parseInt(toks_b[5]);
                                    
                                    if(referenceMap.containsKey(mir_b)){
                                        Integer count_b = referenceMap.get(mir_b);
                                        referenceMap.put(mir_b, count_b+abundance_b);
                                    }else{
                                        referenceMap.put(mir_b, abundance_b);
                                        //First time we have seen a mapping to this mir_id, so increment the number of expressed mirs.
                                        this.expressedReference++;
                                    }
                                }
                                //We know that we have just added the best matches for this sRNA - so any further encounters with this sRNA should be ignored.
                                this.processedSrnaList.put(sRNA, sRNA);
                            }//end what to do if sRNA has best alignments.
                            else{
                                //if nothing was returned, then add this alignmend as it is not in the processed sRNA list and the findBestAlignments returned no best alignment.
                                //So do not add the sRNA to the processed list such to allow all occurances of this sRNA.
                                if(referenceMap.containsKey(mir)){
                                    Integer count = referenceMap.get(mir);
                                    referenceMap.put(mir, count+abundance);
                                }else{
                                    referenceMap.put(mir, abundance);
                                    //First time we have seen a mapping to this mir_id, so increment the number of expressed mirs.
                                    this.expressedReference++;
                                } 
                            }//End else what to do when no best alignments for sRNA found.
                            //We know that the sRNA in this patman alignment is a multimapper so, add the abundance to know that it is a multimapper.
                            multiMapperAbundance.put(sRNA, abundance);
                        }//End if sRNA has previously been seen and 'best matches' processed.          
                    }
                    //Read the next alignment (alignment) of the patman file.
                    line = r.readLine();
                }
            } catch (IOException ex) {
                System.err.println("Error: could not use A provided patman alignment file.");
                Logger.getLogger(SummaryTable.class.getName()).log(Level.SEVERE, null, ex);
            }//end reading the file, all done so store the reference map we have just made.
            this.referenceExpressionMap = referenceMap;
            //SANITY CHECKING! - sum the abundance for the multi-mapping sRNAs.
            Iterator<String> itr = multiMapperAbundance.keySet().iterator();
            while(itr.hasNext()){
                this.multimapper += multiMapperAbundance.get(itr.next());
            }
            //SANITY CHECK! the number of multi-mappers found in this code should match the number of multi-mappers found in the mapping rate code.
            if(this.multimapper != this.totalStranReadsWichHaveTwoOrMoreMappings){
                System.err.println("Error: The total number of multi-mapper counts from two seperate methods do not match (they should!)");
                System.err.println(this.multimapper +" vs "+ this.totalStranReadsWichHaveTwoOrMoreMappings);
            }
            //SANITY CHECK! The number of single matches found in this code should match the number of single matches found in the mapping rate code.
            if(this.mappedOnce != this.totalSrnaReadsWhichHaveExactlyOneMapping){
                System.err.println("Error: The total number of sRNA Reads thave have mapped exactly once does not match the mapping rate code. (It should match!)");
            }
        }
        
        /**
         * Finds the 'best' alignment(s) for the sRNA in the alignment provided and return null if no best alignment(s) were found.
         * Return null logic :- 
         * If there is only a single alignment, then the check should return null because there is no 'best' alignment.
         * If there is multiple alignments all with 0 or all with 1 mismatches, then the check should return null as there is no 'best' alignment.
         * The only time for this method to return an or multiple alignments is if there are a mix of matches with 0 sampleExpressionTable and 1 sampleExpressionTable and the 0 sampleExpressionTable returned.
         * NOTE: Only identified multi-mappers should be passed to this method.
         * @param q_alignment
         * @return Null if there are no better alignments or a list of better alignments.
         * @throws java.io.IOException
         */
        private ArrayList<String> findBestAlignments(String q_alignment) {
            //Get the alignment into its fields.
            String[] q_toks = q_alignment.split("\t");
            String q_miRID = q_toks[0];
            String q_sRNA = q_toks[1].split("\\(")[0];
            int q_abundance = Integer.parseInt( q_toks[1].split("\\(")[1].replace(")", "") );
            int q_mismatches = Integer.parseInt(q_toks[5]);
            
            //Get a list of all alignments for the sRNA in this alignment file.      
            ArrayList<String> alignments = new ArrayList<>();
            //Read the patman alignments file.
            try (BufferedReader r = new BufferedReader(new FileReader(this.patmanAlign))) {
                //Read the first alignment of the patman file.
                String line = r.readLine();
                //While we are not finished reading the patman file.
                while(line != null){
                    //Split the current alignment into its files by a tab delimiter.
                    String[] p_toks = line.split("\t");
                    //Get the sRNA sequence and clean off the bracketed abundance.
                    String p_sRNA = p_toks[1].split("\\(")[0];
                    if(p_sRNA.equalsIgnoreCase(q_sRNA)){
                        alignments.add(line);
                    }
                    line = r.readLine();
                }               
                
            } catch (IOException ex) { 
                System.err.println("Error: can not use a provided patman alignment file.");
                Logger.getLogger(SummaryTable.class.getName()).log(Level.SEVERE, null, ex);
            }
            //Sanity check One - there should not be only one entry in the alignments list as it is a pre-requisite for a multi-mapping alignment to be passed to this method.
            if(alignments.size() <= 1){
                System.err.println("Error: a single matching alignment was passed to the function - findBestAlignments(String q_alignment)");
            }
            ArrayList<String> Zero_MM = new ArrayList<>();
            ArrayList<String> One_MM = new ArrayList<>();
            //Look through the list of all alignments for the sRNA to make a list of best matches.
            for(String align : alignments){
                //Get the alignment into fields.
                String[] c_toks = align.split("\t");
                String c_miRNAID = c_toks[0];
                String c_sRNA = c_toks[1].split("\\(")[0];
                
                //Get the abundance for the patman alignment - we will use this as a sanity check to make 
                //sure the sRNA in each of the alignments is the same abundance (as it should be impossible for a different abundance).
                int c_abundance = Integer.parseInt( c_toks[1].split("\\(")[1].replace(")", "") );
                //Get the number of mismatches for the alignment.
                int c_mismatches = Integer.parseInt(c_toks[5]);              
                
                //Sanity check that each alignment that we look at has the same abundance as the alignment provided.
                if(c_abundance != q_abundance){
                    System.err.println("Error: a sRNA with the same sequence does not hold the same abundance.");
                } 
                //BEST MATCH LOGIC!
                //Look to see if there are multiple alignments and if one or more have 0 mismatches and one or more have 1 mismatches
                switch (c_mismatches) {
                    case 0:
                        Zero_MM.add(align);
                        
                        break;
                    case 1:
                        One_MM.add(align);
                        
                        break;
                    default:
                        System.err.println("Error: an alignment was checked that had more than 1 mismatch in the alignment - this means that more than one mismatch is within the patman alignment record.");
                        break;
                }//end switch                 
            }//end looking through the list.
            if(Zero_MM.size() > 0 && One_MM.isEmpty()){
                //Test to see if the list contained duplicate hits on a single reference.
                ArrayList<String> noDups = this.hasDuplicatesRemoved(Zero_MM);
                if(noDups != null){
                    //It did, so return the list with one hit on a single reference assigned, and any other 0mm hits on other refs for this srna.
                    return noDups;
                } else {
                    //Else if we have a list of only exact matches between different ref seqs, then there is no 'best match' and return null.
                    return null;
                }
            }
            else if(Zero_MM.isEmpty() && One_MM.size() > 0){

                //Test to see if the list contained duplicate hits on a single reference.
                ArrayList<String> noDups = this.hasDuplicatesRemoved(One_MM);
                if(noDups != null){
                    //It did, so return the list with one hit on a single reference assigned, and any other 1mm hits on other refs for this srna.
                    return noDups;
                } else {
                    //If we have a list of only One Mismatch accross multiple references, then there is no 'best match' and return null.
                    return null;
                }
            }
            else if(!Zero_MM.isEmpty() && !One_MM.isEmpty()){
                //If there are alignments for both 0mm and 1mm, then there are some 'best matches' and return the best (0mm)
                //Either with dups removes (if found) or as listed because there are no dupes.                
                //Test to see if the list contained duplicate hits on a single reference.
                ArrayList<String> noDups = this.hasDuplicatesRemoved(Zero_MM);
                if(noDups != null){
                    //It did, so return the list with one hit on a single reference assigned, and any other 0mm hits on other refs for this srna.
                    return noDups;
                } else {
                    return Zero_MM;
                }
            }
            //We have a sanity check to make sure that there is no other option - this should never happen so report an error!
            System.err.println("Error: an error has occured that is not well defined. A best match may not have been reported.");
            return null; //This return of null should never be used and has been reported as an error.   
        }//end method.
         
        
        public ArrayList<String> hasDuplicatesRemoved(ArrayList<String> alignments){
                //IF we have a list of alignments, if any are on a single reference, we only want to matches on a single reference once,
                //so filter out duplicate mappings of an sRNA single reference sequence,
                //but also retain alignments for that sRNA which map with 0 mm to other reference sequences.
                boolean hasDups = false;
                //Test to see if best alignments list contain multimappers on same reference sequence.
                HashMap<String,String> nonDuplicatedAssignments = new HashMap<>();
                for(String align : alignments){
                    //Split the current best alignment into fields by tab.
                    String[] toks_a = align.split("\t");
                    //Get the current best mir_id.
                    String mir_a = toks_a[0];
                    //We already know that all sRNA seqs in the alignments list are the same.
                    //So we only deal with any duplicate reference IDs by using collisions within the map to overwrite duplication
                    //therefore assigning an arbitrary position alignment (last-in).
                    if(nonDuplicatedAssignments.containsKey(mir_a)){
                        hasDups = true;
                    }
                    nonDuplicatedAssignments.put(mir_a, align);   
                }
                
                if(hasDups){
                    //Then we found duplicate 0mm alignments on a single ref and we should assign one of them as an arbitrary best match.
                    //We have a map of alignments which are non-duplicated for a single reference seq, so make a list
                    //and continue allocation of the best (arbitrary) alignment.
                    ArrayList<String> Zero_MM_noDup = new ArrayList<>();
                    Iterator<String> itr = nonDuplicatedAssignments.keySet().iterator();
                    while(itr.hasNext()){
                        String refID = itr.next();
                        String alignment = nonDuplicatedAssignments.get(refID);
                        Zero_MM_noDup.add(alignment);
                    }
                    return Zero_MM_noDup;
                }
                return null;
        }
        
        /**
         * Method to check if a query sRNA is a multi-mapper.
         * NOTE: If a value of 1 is received by a calling method, then the calling method must
         * decide how to deal with this sRNA - at best this method tells the calling method that
         * the sRNA provided has multiple mappings to one or more miRNA.
         * For example, the calling method should test if the alignments for this sRNA contains 0 mismatches. 
         * If there are, then only the 0 mismatch alignments should be allocated and mappings with 1 mismatch should be discarded. 
         * However, if there are only alignments with 1 mismatch then all alignments should be retained.
         * @param sRNAToTest - The sRNA query.
         * @return 0 = exactly 1 hit on 1 miRNA. 1 = sRNA is a multi-mapper.
         * @throws FileNotFoundException
         * @throws IOException 
         */
        public int checkMultimapperStatus(String sRNAToTest) throws FileNotFoundException, IOException{
            //Assigned value for the mapping status. An assigned value of 99, as 99 will be an error code where no status was assigned.
            int assignedVal = 99; 
            //Hashmap <miRNA:PatmanHitCount> to hold the number of patman algnments have been made for the miRNAs.
            //Note that because the sRNAs are in NR format, we are storing only the number of times an alignment
            //was found, and not the number of reads that were aligned.
            HashMap<String, Integer> hitsCount;
            //A list to hold the alignments (lines) read in from the patman file.
            ArrayList<String> alignments;
            //Read the patman alignments file.
            try (BufferedReader r = new BufferedReader(new FileReader(this.patmanAlign))) {
                //Initialize hash and list.
                hitsCount = new HashMap<>();
                alignments = new ArrayList<>();
                //Read the first alignment of the patman file.
                String line = r.readLine();
                //While we are not finished reading the patman file.
                while(line != null){
                    //Split the current alignment into its files by a tab delimiter.
                    String[] toks = line.split("\t");
                    //Get the miRNA id which this alignment (where a sRNA has mapped to a miRNA).
                    String mir = toks[0];
                    //Get the sRNA sequence and clean off the bracketed abundance.
                    String sRNA = toks[1].split("\\(")[0];
                    //IMPORTANT: We are only going to do anything for an alignment if it matches the sRNA being queried. So,
                    //If the sRNA sequence currently being queried is exactly the same as the current sRNA being read in the patman file 
                    if(sRNA.equalsIgnoreCase(sRNAToTest)){
                        //Then add the current alignment to the list of alignments (patman file lines) to keep track of alignments for this sRNA sequence.
                        alignments.add(line);
                        //If we have seen an alignment hit in the patman file for this miR before, sampleIndex.e. the query sRNA maps multiple times to this mir_id.
                        if(hitsCount.containsKey(mir)){
                            //Then get the number of times this sRNA has aligned to this miRNA before.
                            int numHits = hitsCount.get(mir);
                            //And increment the number of times this sRNA has aligned to this miRNA by 1 more.
                            hitsCount.put(mir, numHits + 1);
                        }else{
                            //Else, this is the first time we have had a hit for this sRNA on this miRNA, so add the mir_id and its first count.
                            hitsCount.put(mir, 1);
                        }
                    }
                    //Move on to the next alignment in the patman file.
                    line = r.readLine();
                }//end while.
            }
            //All done looking for alignments for the query sRNA, so now to look at what the mapping status of the query sRNA.
            
            //If the number of ref miRNAs that we had a hit for the query sRNAs was only one, and the number of hits for the query sRNA on that
            //miRNA was only 1.
            if(hitsCount.size() == 1 && alignments.size() == 1 ){
                //Then we have found Exactly one hit on only one miRNA for the query sRNA.
                assignedVal = 0;
            }
            //Else if we have only 1 ref miRNA and the number of hits on that miRNA is equal or more than 2 OR
            //we have multiple miRNAs in the list.
            else if( (hitsCount.size() == 1 && alignments.size() >= 2) || hitsCount.size() > 1){
                //Then the query sRNA is a multi-mapper and the best alignments should be identified.
                assignedVal = 1;
            }
            //If no assignment of single exact match or multi-mapper was assigned.
            if(assignedVal == 99){
                //Then some error has occured - this should never be possible.
                System.err.println("Error: An sRNA was not assigned a mapping status.");
            }
            return assignedVal;
        }//end method.
        
        /**
         * Calculates the mapping rates for the sRNAs in this sample.
         */
        private void calculateMappingRate() {
            
            //Local variable to hold the mapping counts.
            Integer local_totalSrnasWhichHaveMapping = 0;
            double local_totalSrnasWhichHaveMapping_PERCENT;      
            Integer local_totalSrnasWhichHaveExactlyOneMapping = 0;
            double local_totalSrnasWhichHaveExactlyOneMapping_PERCENT; 
            Integer local_totalStransWichHaveTwoOrMoreMappings = 0;
            double local_totalStransWichHaveTwoOrMoreMappings_PERCENT;  
            Integer local_totalSrnasWithoutMapping = 0;
            double local_totalSrnasWithoutMapping_PERCENT;
            //A hashmap containing Srna Mapping Locations counts <sRNA_sequence:mapping counts>.
            HashMap<String, Integer> sRNA_Mapping_Location_Count = new HashMap<>();
            //A hashmap containing sRNA read counts <sRNA_sequence:read counts>.
            HashMap<String, Integer> sRNA_Read_Counts = new HashMap<>();
            
            //Get sRNA Read Counts in to a hashmap and initialize a hashmap of mapping locations for each sRNA (initalized to zero locations).
            //Open the sRNA fasta file for this sample.
            try (BufferedReader r = new BufferedReader(new FileReader(this.SRNA_FASTA))) {
                //Read the first two lines - header and sequence lines.
                String FastaHeaderLine = r.readLine();
                String FastaSequenceLine = r.readLine(); //(sequence alignment in fasta is never used).
                //While there are lines to read.
                while(FastaHeaderLine != null){
                    //Get only the sequence from the fasta header alignment - (clean off the abundance and > symbol)
                    String srnaSequence_key = FastaHeaderLine.split("\\(")[0].replace(">", "");
                    //Get the abundance for the sRNA sequence.
                    Integer abundance = Integer.parseInt(FastaHeaderLine.split("\\(")[1].replace(")", "") );
                    //SANITY CHECK - make sure that this fasta file is NR as expected, i.e. >[UNIQUESrnaSequence](abundance).
                    if(sRNA_Mapping_Location_Count.containsKey(srnaSequence_key)){
                        System.err.println("Input Error: There was a collision in a hashmap for two sequences of the same composition. This Implies that the provided sRNA file is not NR.");
                        System.exit(1);
                    } else {
                        //This else block should ALWAYS be entered into.
                        //So add the sRNA sequence as a key to the mapping location count hashmap.
                        sRNA_Mapping_Location_Count.put(srnaSequence_key, 0);
                        //Also add the same key and its abundance to the sRNA read counts hashmap.
                        sRNA_Read_Counts.put(srnaSequence_key, abundance);
                    }
                    //Read the next entry in the fasta file.
                    FastaHeaderLine = r.readLine();
                    FastaSequenceLine = r.readLine();
                }
            } catch (IOException ex) {
                System.err.println("Error: could not use provided srna fasta file.");
                Logger.getLogger(SummaryTable.class.getName()).log(Level.SEVERE, null, ex);
            }//End reading the fasta file.
             
            //Find the number of times each sRNA sequence has a mapping to the reference sequence and keep a record of it in sRNA_Mapping_Location_Count
            //i.e. If an sRNA has a mapping, then add it to the mapping count for the srna in the mapping counts hashmap.
            //Open a stream to the patman file.
            try (BufferedReader r = new BufferedReader(new FileReader(this.getPatmanAlign()))) {
                //Read the first alignment.
                String alignment = r.readLine();
                //While we have lines to read.
                while(alignment != null){
                    //Split the alignment record up into fields on tabs.
                    String[] toks = alignment.split("\t");
                    //Extract the sequence and clean off the abundance and > symbol.
                    String sequence = toks[1].split("\\(")[0].replace(">", "");
                    //If the sRNA has already been seen i.e. it already has a mapping.
                    if(sRNA_Mapping_Location_Count.containsKey(sequence)){
                        //Then get the existing number of mappings found for the sRNA.
                        Integer count = sRNA_Mapping_Location_Count.get(sequence);
                        //And increment it by 1 as another mapping has been found.
                        sRNA_Mapping_Location_Count.put(sequence, count+1);
                    }else{
                        //Else, this is the first mapping we have seen for this sRNA, so record it in the hashmap for this sRNA.
                        sRNA_Mapping_Location_Count.put(sequence, 1);
                    }
                    //Read the next patman alignment.
                    alignment = r.readLine();
                }//end while.
            } catch (IOException ex) {
                System.err.println("Error: could not use provided patman alignment file.");
                Logger.getLogger(SummaryTable.class.getName()).log(Level.SEVERE, null, ex);
            }//end the reading of the patman file.
            
            //Get the grand sum of the mapping counts for sRNAs which either A) did not map at all, B) mapped exactly once, and C) mapped more than once.
            //Get the key set for the sRNA mapping location count hashmap - the key is the sRNA sequence (NOTE: each and every sRNA is in the hashmap).
            Iterator<String> itr = sRNA_Mapping_Location_Count.keySet().iterator();
            //While we have sRNAs in the hashmap to be processed.
            while(itr.hasNext()){
                //Get the next sRNA key.
                String key = itr.next();
                //Get the number of times the sRNA has had a mapping.
                Integer count = sRNA_Mapping_Location_Count.get(key);
                //If the sRNA has no mapping at all.
                if(count == 0){
                    //Then increment the local variable holding NO MAPPING by the number of reads for the sRNA sequence.
                    local_totalSrnasWithoutMapping += sRNA_Read_Counts.get(key);
                }
                //If the sRNA has exactly one mapping.
                else if(count == 1){
                    //Then increment the local variable holding exactly one mapping by the number of reads for the sRNA sequence.
                    local_totalSrnasWhichHaveExactlyOneMapping += sRNA_Read_Counts.get(key);
                    //Also, as we have a mapping, we increment the local variable which keeps count of all reads which have a mapping by the number of reads for the sRNA sequence.
                    local_totalSrnasWhichHaveMapping += sRNA_Read_Counts.get(key);
                }
                //If the sRNA has two ore more mappings.
                else if(count >= 2){
                    //Then we increment thelocal variable holding counts of sRNA which have two ore more mappings by the number of reads for the sRNA sequence.
                    local_totalStransWichHaveTwoOrMoreMappings += sRNA_Read_Counts.get(key);
                    //Also, as we have a mapping, we increment the local variable which keeps count of all reads which have a mapping by the number of reads for the sRNA sequence.
                    local_totalSrnasWhichHaveMapping += sRNA_Read_Counts.get(key);
                }
                //Else, SANITY CHECK! something has gone wrong because there should never be a sRNA which has not been counted as either not mapping, mapping once or mapping twice or more.
                else{
                    //So tell somebody that something weird has happend and what the effect was.
                    System.err.println("Warning: A srna mapping has been found which remains unplaced and has been discarded!");
                } 
            }
            
            //Based on the counts we have just obtained - calculate the counts as a percent of total sRNA reads in the fasta file.
            local_totalSrnasWhichHaveMapping_PERCENT = ((double)local_totalSrnasWhichHaveMapping / (double)this.totalSrnasInSample) * (double)100.00;
            local_totalSrnasWhichHaveExactlyOneMapping_PERCENT = ((double)local_totalSrnasWhichHaveExactlyOneMapping / (double)this.totalSrnasInSample) * (double)100.00;
            local_totalStransWichHaveTwoOrMoreMappings_PERCENT = ((double)local_totalStransWichHaveTwoOrMoreMappings / (double)this.totalSrnasInSample) * (double)100.00;
            local_totalSrnasWithoutMapping_PERCENT = ((double)local_totalSrnasWithoutMapping / (double)this.totalSrnasInSample) * (double)100.00;

            //Place all of the local fields we have been working with into a class field.
            this.totalSrnasReadsWhichHaveMapping = local_totalSrnasWhichHaveMapping;
            this.totalSrnaReadsWhichHaveMapping_PERCENT = local_totalSrnasWhichHaveMapping_PERCENT; 
            this.totalSrnaReadsWhichHaveExactlyOneMapping = local_totalSrnasWhichHaveExactlyOneMapping;
            this.totalSrnaReadsWhichHaveExactlyOneMapping_PERCENT = local_totalSrnasWhichHaveExactlyOneMapping_PERCENT;    
            this.totalStranReadsWichHaveTwoOrMoreMappings = local_totalStransWichHaveTwoOrMoreMappings;
            this.totalStranReadsWichHaveTwoOrMoreMappings_PERCENT = local_totalStransWichHaveTwoOrMoreMappings_PERCENT;
            this.totalSrnaReadsWithoutMapping = local_totalSrnasWithoutMapping;
            this.totalSrnaReadsWithoutMapping_PERCENT = local_totalSrnasWithoutMapping_PERCENT;  
        
        }

        public String getSampleName() {
            return sampleName;
        }

        public File getSRNA_FASTA() {
            return SRNA_FASTA;
        }

        public Integer getTotalSrnasInSample() {
            return totalSrnasInSample;
        }

        public File getPatmanAlign() {
            return patmanAlign;
        }

        public Integer getTotalSrnasWhichHaveMapping() {
            return totalSrnasReadsWhichHaveMapping;
        }

        public Integer getTotalSrnasWhichHaveExactlyOneMapping() {
            return totalSrnaReadsWhichHaveExactlyOneMapping;
        }

        public Integer getTotalStransWichHaveTwoOrMoreMappings() {
            return totalStranReadsWichHaveTwoOrMoreMappings;
        }

        public Integer getTotalSrnasThatDidNotMap() {
            return totalSrnaReadsWithoutMapping;
        }

        public double getTotalSrnasWhichHaveMapping_PERCENT() {
            return totalSrnaReadsWhichHaveMapping_PERCENT;
        }

        public double getTotalSrnasWhichHaveExactlyOneMapping_PERCENT() {
            return totalSrnaReadsWhichHaveExactlyOneMapping_PERCENT;
        }

        public double getTotalStransWichHaveTwoOrMoreMappings_PERCENT() {
            return totalStranReadsWichHaveTwoOrMoreMappings_PERCENT;
        }

        public double getTotalSrnasWithoutMapping_PERCENT() {
            return totalSrnaReadsWithoutMapping_PERCENT;
        }

        public String getCondition() {
            return condition;
        }

        public String getReplicate() {
            return replicate;
        }
        
        public int getExpressedReferenceCount() {
            return expressedReference;
        }

        private HashMap<String, Integer> getSampeReferenceExpressionTable() {
            return this.referenceExpressionMap;
        }
    
    }//end inner class.
      
}//end outer class.
