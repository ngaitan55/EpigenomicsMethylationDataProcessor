package HelperScripts;
import java.util.Arrays;
import java.util.List;

import ngsep.genome.ReferenceGenome;
import ngsep.sequences.QualifiedSequence;
import ngsep.sequences.QualifiedSequenceList;

public class FastaCountMethilationContexts {

	public static void main(String[] args) throws Exception {
		String referenceFile = args[0];
		int windowLength = Integer.parseInt(args[1]);
		ReferenceGenome genome = new ReferenceGenome(referenceFile);
		QualifiedSequenceList sequences = genome.getSequencesList();
		for(QualifiedSequence qseq:sequences) {
			String seqChars = qseq.getCharacters().toString();
			int n = seqChars.length();
			for(int i=0;i<n;i+=windowLength) {
				int [] contextCounts = new int [3];
				Arrays.fill(contextCounts, 0);
				for(int j=0;j<windowLength;j++) {
					int next = i+j;
					if(next>=n-2) break;
					char bp0 = seqChars.charAt(next);
					if(bp0=='C') {
						char bp1 = seqChars.charAt(next+1);
						char bp2 = seqChars.charAt(next+2);
						if(bp1=='G') contextCounts[0]++;
						else if (bp1!='C') {
							if(bp2=='G') contextCounts[1]++;
							else if (bp2!='C') contextCounts[2]++;
						}
					} else if (bp0=='G' && next>1) {
						//Negative strand
						char bp1 = seqChars.charAt(next-1);
						char bp2 = seqChars.charAt(next-2);
						if(bp1=='C') contextCounts[0]++;
						else if (bp1!='G') {
							if(bp2=='C') contextCounts[1]++;
							else if (bp2!='G') contextCounts[2]++;
						}	
					}		
				}
				System.out.println(""+qseq.getName()+"\t"+(i+1)+"\t"+(i+windowLength)+"\t"+contextCounts[0]+"\t"+contextCounts[1]+"\t"+contextCounts[2]);
			}
		}
	}
}
