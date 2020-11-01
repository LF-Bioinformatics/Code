package makefastanr;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

public class MakeFastaNR {

    private final File inFile;
    private HashMap<String, String> nrSequence = new HashMap<>();
    
    public MakeFastaNR(File in){
        this.inFile = in;
    }
    
    public void makeNR() throws FileNotFoundException, IOException{
        BufferedReader r = new BufferedReader(new FileReader(inFile));
        String lineHead = r.readLine();
        String lineSeq = r.readLine();
        while(lineSeq != null){
            if( this.nrSequence.containsKey(lineSeq) ){
                String newHead = this.nrSequence.get(lineSeq) + "_" + lineHead;
                this.nrSequence.put(lineSeq, newHead);
            }else{
                this.nrSequence.put(lineSeq, lineHead);
            }
            lineHead = r.readLine();
            lineSeq = r.readLine();
        }
    }
    
    public void condenseSubsequences(){  
        HashMap<String, String> nrSubSequence = new HashMap<>();
        HashMap<String, String> nrFinal = new HashMap<>();
        nrSubSequence.putAll(this.nrSequence);
        nrFinal.putAll(this.nrSequence);
        Iterator<String> itr = this.nrSequence.keySet().iterator();
        while(itr.hasNext()){
            String mirSeq_1 = itr.next();
            String mirID_1 = this.nrSequence.get(mirSeq_1);
            Iterator<String> local = nrSubSequence.keySet().iterator();
            while(local.hasNext()){
                String mirSeq_2 = local.next();
                String mirID_2 = nrSubSequence.get(mirSeq_2);
                if(mirSeq_1.contains(mirID_2)){
                    //append mir 1 id to show that mir 1 also contains mir 2 sequence and remove mir 2
                    String newID = mirID_1 + " " + mirID_2;
                    nrSubSequence.put(mirSeq_1, newID);
                    nrSubSequence.remove(mirSeq_2);
                }else if(mirSeq_2.contains(mirID_1)){
                    //append mir 2 id to show that mir 2 also contains mir 1 sequence and remove mir 1
                    String newID = mirID_2 + " " + mirID_1;
                    nrSubSequence.put(mirSeq_2, newID);
                    nrSubSequence.remove(mirSeq_1);
                }  
            }
        }
        this.nrSequence = nrSubSequence;
        
    }
    
    public void outputNR(){
        Iterator<String> itr = this.nrSequence.keySet().iterator();
        while(itr.hasNext()){
            String seq = itr.next();
            String head = this.nrSequence.get(seq);
            System.out.println(head);
            System.out.println(seq);
        }
    }
    
    public static void main(String[] args) throws IOException {
        
        if(args.length == 0){
            System.out.println("Java -Jar MakeFastaNR.jar [path to input fasta file]");
            System.out.println("Input fasta file will be made NR and output to STD out.");
        }else{
            MakeFastaNR nr = new MakeFastaNR(new File(args[0]));
            nr.makeNR();
            nr.condenseSubsequences();
            nr.outputNR();
        }
    }    
}