package org.geatools.data.structure;

import java.util.Comparator;

public class SeqEncode {
   public int seqID;
   public long numEncode;
   public long numRevEncode;

   public boolean isDup=false;
   public int dupNum=-1;
   
   public static class CompSeqEncode implements Comparator<SeqEncode> {
			
		private int mod = 1;
		public CompSeqEncode(boolean desc) {
		  if (desc) mod =-1;
		}
			        
		public int compare(SeqEncode  seq1, SeqEncode  seq2){
		  return  mod*Long.valueOf(seq1.numEncode).compareTo(seq2.numEncode);
		}
			
	}
   
	public static class CompSeqID implements Comparator<SeqEncode> {
			
		private int mod = 1;
		public CompSeqID(boolean desc) {
		  if (desc) mod =-1;
		}
			        
		public int compare(SeqEncode  seq1, SeqEncode  seq2){
		  return  mod*Integer.valueOf(seq1.seqID).compareTo(seq2.seqID);
		}
			
	 }

}
