package expressionmatrix;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Generate expression matrix from PatMaN alignments.
 */
public class ExpressionMatrix {

    private String outputPrefix;
    private File sampleTable;
    private SummaryTable summaryTable;
    private String[] sampleName;
    private String[] condition;
    private File[] fastaNonRedundant;	
    private File[] fastq;	
    private String[] replicateName;
    private File[] patmanAlignment;

    /**
     * Constructs an empty Expression Matrix
     */
    public ExpressionMatrix() { }
    
    /**
     * Loads the sample table information into the expression matrix.
     * @param args positional arguments [0]=sample table path, [1]=output path prefix.
     */
    public void loadSampleTable(String[] args) {
        this.sampleTable = new File(args[0]);
        this.outputPrefix = args[1];
        ArrayList<String> tempSampleTable = new ArrayList<>();
        try (BufferedReader r = new BufferedReader(new FileReader(sampleTable))) {
            r.readLine();
            String line = r.readLine();
            while(line != null){ tempSampleTable.add(line); line = r.readLine(); }
        } catch (IOException ex) {
            Logger.getLogger(ExpressionMatrix.class.getName()).log(Level.SEVERE, null, ex);
        }
        this.sampleName = new String[tempSampleTable.size()];
        this.condition = new String[tempSampleTable.size()];
        this.replicateName = new String[tempSampleTable.size()];
        this.fastaNonRedundant = new File[tempSampleTable.size()];
        this.fastq = new File[tempSampleTable.size()];
        this.patmanAlignment = new File[tempSampleTable.size()];
        for(int i = 0; i < tempSampleTable.size(); i++){
            String[] toks = tempSampleTable.get(i).split("\t");
            this.sampleName[i] = toks[0];
            this.condition[i] = toks[1];
            this.replicateName[i] = toks[2];
            if(toks[4].endsWith("/")){ this.fastaNonRedundant[i] = new File(toks[4]+toks[3]); }
            else{ this.fastaNonRedundant[i] = new File(toks[4]+"/"+toks[3]); }
            if(toks[6].endsWith("/")){ this.fastq[i] = new File(toks[6]+toks[5]); }
            else{ this.fastq[i] = new File(toks[6]+"/"+toks[5]); }
            if(toks[8].endsWith("/")){ this.patmanAlignment[i] = new File(toks[8]+toks[7]); }
            else{ this.patmanAlignment[i] = new File(toks[8]+"/"+toks[7]); }
        }
        System.out.println("Number of samples in experiment:\t"+tempSampleTable.size());
    }
    
    /**
     * Generates and expression matrix and summary table document.
     */
    public void generateExpressionMatrixSummary() {
        this.summaryTable = new SummaryTable(this.fastaNonRedundant, this.patmanAlignment, 
                this.sampleName, this.condition, this.replicateName, this.outputPrefix);
        this.summaryTable.GenerateExperimentSummary();
    }
    
    /**
     * Launch the program with the provided command line arguments. 
     * @param args positional arguments required [0]=sample table in tsv format, [1] output prefix. 
     */
    public static void main(String[] args) {
        
        if(args == null || args.length == 0){
            System.out.println("ExpressionMatrix version 1.2");
            System.out.println("Usage: java -Jar ExpressionMatrix.jar [sample table in .tsv format] [output prefix]");
            System.exit(0);           
        }
        ExpressionMatrix em = new ExpressionMatrix();
        em.loadSampleTable(args);
        em.generateExpressionMatrixSummary();
    }    
}
