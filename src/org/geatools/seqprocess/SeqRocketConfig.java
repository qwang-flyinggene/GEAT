package org.geatools.seqprocess;

import java.io.File;
import java.util.List;

import org.geatools.data.structure.SeqCompoRecognizer;
import org.geatools.data.structure.SeqRocket;
import org.geatools.operation.SeqOperation;

public class SeqRocketConfig extends SeqRocketRecognition {
	
	int barExactMaxStart=2;	
	double barMinAlignRatio=0.8d;
	double barMaxMismatchRatio=0.1d;
	double barMaxGapRatio=0.1d;	  
	  
	double primerMinAlignRatio=0.7d;
	double primerMaxMismatchRatio=0.1d;
	double primerMaxGapRatio=0.1d;
	int primerLeftShift=5;
	  
	double primerContMinAlignRatio=0.8d;
	double primerContMaxMismatchRatio=0.1d;
	double primerContMaxGapRatio=0.1d;
	  
	double baitMinAlignRatio=0.25d;
	double baitMaxMismatchRatio=0.2d;
	double baitMaxGapRatio=0.2d;
	//int baitLeftShift=10;
	int baitLeftShift=0;
	  
	int baitBrkMinAlign=12;
	double baitBrkMinAlignRatio=0.5d;
	double baitBrkMaxMismatchRatio=0.2d;
	double baitBrkMaxGapRatio=0.2d;
	//int baitBrkLeftShift=10;
	int baitBrkLeftShift=0;
	  
	int baitArmMinAlign=20;
	double baitArmMinAlignRatio=0.6d;
	double baitArmMaxMismatchRatio=0.2d;
	double baitArmMaxGapRatio=0.2d;
	//int baitArmLeftShift=10;  
	int baitArmLeftShift=0; 
	  
	int rc3MinAlignLen=16;
	double rc3MinAlignRatio=0.2d;	
	double rc3MaxMismatchRatio=0.2d;
	double rc3MaxGapRatio=0.1d;
	double rc3MaxTStartRatio=0.75d;
	//int rc3MaxQStar=10;
	//int rc3RightShift=0;
	  
	double rc3PrimerContMinAlignRatio=0.3d;
	double rc3PrimerContMaxMismatchRatio=0.2d;
	double rc3PrimerContMaxGapRatio=0.2d;
	int rc3PrimerContRightShift=0;
	  
	double rc3PrimerMinAlignRatio=0.3d;
	double rc3PrimerMaxMismatchRatio=0.2d;
	double rc3PrimerMaxGapRatio=0.2d;
	int rc3PrimerRightShift=0;
	  
	double rc3BarcodeMinAlignRatio=0.5d;
	double rc3BarcodeMaxMismatchRatio=0.2d;
	double rc3BarcodeMaxGapRatio=0.2d;
	int rc3BarcodeRightShift=0;
	  
	double rc3SeqLinkerMinAlignRatio=0.3d;
	double rc3SeqLinkerMaxMismatchRatio=0.2d;
	double rc3SeqLinkerMaxGapRatio=0.2d;
	int rc3SeqLinkerRightShift=0;	  			 

	boolean mask_left_barcode=true;
	boolean mask_left_primer=true;
	boolean mask_left_primerCont=true;
	boolean mask_left_bait=true;
	boolean mask_left_baitBrk=true;
	boolean mask_left_baitArm=true;   
	boolean mask_right_rc3together=true;
	boolean mask_right_rc3primerCont=true;
	boolean mask_right_rc3primer=true;
	boolean mask_right_rc3barcode=true;
	boolean mask_right_rc3seqLinker=true;
	  
	boolean trim_left_barcode=true;
	boolean trim_left_primer=true;
	boolean trim_left_primerCont=true;
	boolean trim_left_bait=true;
	boolean trim_left_baitBrk=true;
	boolean trim_left_baitArm=true; 
	boolean trim_right_rc3together=true;
	boolean trim_right_rc3primerCont=true;
	boolean trim_right_rc3primer=true;
	boolean trim_right_rc3barcode=true;
	boolean trim_right_rc3seqLinker=true;	  

	boolean saveFinalSeq=true;
	boolean saveFinalSeqAsFASTA=true;
	boolean saveFinalSeqAsFASTQ=false;
	boolean saveFinalSeqAsHTML=true;
	  /*
	boolean barcode_save=true;
	boolean primer_save=false;
	boolean primerCont_save=false;
	boolean bait_save=false;
	boolean baitBrk_save=false;
	boolean baitArm_save=false;
	  */
	boolean barcode_saveAsFASTA=false;
	boolean barcode_saveAsFASTQ=false;
	boolean primer_saveAsFASTA=false;
	boolean primer_saveAsFASTQ=false;
	boolean primerCont_saveAsFASTA=false;
	boolean primerCont_saveAsFASTQ=false;
	boolean bait_saveAsFASTA=false;
	boolean bait_saveAsFASTQ=false;
	boolean baitBrk_saveAsFASTA=false;
	boolean baitBrk_saveAsFASTQ=false;
	boolean baitArm_saveAsFASTA=false;
	boolean baitArm_saveAsFASTQ=false;	
	  
	public SeqRocketConfig(){ 

	}
		  
	public void setBarcodeMaxExactStart(int start){ 
	   if(start<1) start=1;
	   barExactMaxStart=start;
	} 
	  
	public int getBarcodeMaxExactStart(){ 
	   return barExactMaxStart;
	} 
	  
	public void setRC3MaxTStartPCTL(double ratio){ 
	   rc3MaxTStartRatio=ratio;
	} 

    public double getRC3MaxTStartPCTL(){ 
	   return rc3MaxTStartRatio;
	}
	  
	public void setRC3MinAlignLen(int minLen){ 
	   rc3MinAlignLen=minLen;
	} 

	public int getRC3MinAlignLen(){ 
	   return rc3MinAlignLen;
	}
	  
	public void setMaskLeft(String name, boolean isDoMask){ 
		  
		  mask_left_barcode=false;
		  mask_left_primer=false;
		  mask_left_primerCont=false;
		  mask_left_bait=false;
		  mask_left_baitBrk=false;
		  mask_left_baitArm=false;		
		  
		  if(name==null) {
			  mask_left_barcode=isDoMask;
			  mask_left_primer=isDoMask;
			  mask_left_primerCont=isDoMask;
			  mask_left_bait=isDoMask;
			  mask_left_baitBrk=isDoMask;
			  mask_left_baitArm=isDoMask;
		  }else if(name.equalsIgnoreCase(BARCODE_NAME_DEFINITION)) {		
			  mask_left_barcode=isDoMask; 
		  }else if(name.equalsIgnoreCase(PRIMER_NAME_DEFINITION)) {
			  mask_left_primer=isDoMask;	 
		  }else if(name.equalsIgnoreCase(PRIMERCONT_NAME_DEFINITION)) {	  		   
			  mask_left_primerCont=isDoMask;	 
		  }else if(name.equalsIgnoreCase(BAIT_NAME_DEFINITION)) {	  			
			  mask_left_bait=isDoMask;		    
		  }else if(name.equalsIgnoreCase(BAITBRK_NAME_DEFINITION)) {	  			
			  mask_left_baitBrk=isDoMask;			  
		  }else if(name.equalsIgnoreCase(BAITARM_NAME_DEFINITION)) {
			  mask_left_baitArm=isDoMask;		 
		  }

	}
	  
	public void setTrimLeft(String name, boolean isDoTrim){ 
		  
		  trim_left_barcode=false;
		  trim_left_primer=false;
		  trim_left_primerCont=false;
		  trim_left_bait=false;
		  trim_left_baitBrk=false;
		  trim_left_baitArm=false;		
		  
		  if(name==null) {
			  trim_left_barcode=isDoTrim;
			  trim_left_primer=isDoTrim;
			  trim_left_primerCont=isDoTrim;
			  trim_left_bait=isDoTrim;
			  trim_left_baitBrk=isDoTrim;
			  trim_left_baitArm=isDoTrim;		  
		  }else if(name.equalsIgnoreCase(BARCODE_NAME_DEFINITION)) {		
			  trim_left_barcode=isDoTrim; 
		  }else if(name.equalsIgnoreCase(PRIMER_NAME_DEFINITION)) {
			  trim_left_primer=isDoTrim;	 
		  }else if(name.equalsIgnoreCase(PRIMERCONT_NAME_DEFINITION)) {	  		   
			  trim_left_primerCont=isDoTrim;	 
		  }else if(name.equalsIgnoreCase(BAIT_NAME_DEFINITION)) {	  			
			  trim_left_bait=isDoTrim;		    
		  }else if(name.equalsIgnoreCase(BAITBRK_NAME_DEFINITION)) {	  			
			  trim_left_baitBrk=isDoTrim;			  
		  }else if(name.equalsIgnoreCase(BAITARM_NAME_DEFINITION)) {
			  trim_left_baitArm=isDoTrim;		 
		  }

	}

	  
	public void setMaskRight(String name, boolean isDoMask){ 
		   
		  mask_right_rc3together=false;
		  mask_right_rc3primerCont=false;	  
		  mask_right_rc3primer=false;
		  mask_right_rc3barcode=false;
		  mask_right_rc3seqLinker=false;	
		  if(name==null) {
			  mask_right_rc3together=isDoMask;
			  mask_right_rc3primerCont=isDoMask;	  
			  mask_right_rc3primer=isDoMask;
			  mask_right_rc3barcode=isDoMask;
			  mask_right_rc3seqLinker=isDoMask;	 
		  }else if(name.equalsIgnoreCase(RC3TOGETHER_NAME_DEFINITION)) {		
			  mask_right_rc3together=isDoMask; 
		  }else if(name.equalsIgnoreCase(RC3BARCODE_NAME_DEFINITION)) {		
			  mask_right_rc3barcode=isDoMask; 
		  }else if(name.equalsIgnoreCase(RC3PRIMER_NAME_DEFINITION)) {
			  mask_right_rc3primer=isDoMask;	 
		  }else if(name.equalsIgnoreCase(RC3PRIMERCONT_NAME_DEFINITION)) {	  		   
			  mask_right_rc3primerCont=isDoMask;	 
		  }else if(name.equalsIgnoreCase(RC3SEQLINKER_NAME_DEFINITION)) {	  		   
			  mask_right_rc3primerCont=isDoMask;	 
		  }

	}
	 
	public void setTrimRight(String name, boolean isDoTrim){ 
		   
		  trim_right_rc3together=false;
		  trim_right_rc3primerCont=false;	  
		  trim_right_rc3primer=false;
		  trim_right_rc3barcode=false;
		  trim_right_rc3seqLinker=false;	
		  
		  if(name==null) {
			 trim_right_rc3together=isDoTrim;
			 trim_right_rc3primerCont=isDoTrim;	  
			 trim_right_rc3primer=isDoTrim;
			 trim_right_rc3barcode=isDoTrim;
			 trim_right_rc3seqLinker=isDoTrim;	
		  }else if(name.equalsIgnoreCase(RC3TOGETHER_NAME_DEFINITION)) {		
			 trim_right_rc3together=isDoTrim; 
		  }else if(name.equalsIgnoreCase(RC3BARCODE_NAME_DEFINITION)) {		
			 trim_right_rc3barcode=isDoTrim; 
		  }else if(name.equalsIgnoreCase(RC3PRIMER_NAME_DEFINITION)) {
			 trim_right_rc3primer=isDoTrim;	 
		  }else if(name.equalsIgnoreCase(RC3PRIMERCONT_NAME_DEFINITION)) {	  		   
			 trim_right_rc3primerCont=isDoTrim;	 
		  }else if(name.equalsIgnoreCase(RC3SEQLINKER_NAME_DEFINITION)) {	  		   
			 trim_right_rc3primerCont=isDoTrim;	 
		  }

	}	  

	public void setSaveRecognizedSeqFormat(boolean saveAsFASTA, boolean saveAsFASTQ){ 	  
		
		  saveFinalSeq=false;
		  saveFinalSeqAsFASTA=false;
		  saveFinalSeqAsFASTQ=false;
		  if(saveAsFASTA || saveAsFASTQ) {
			 saveFinalSeq=true;
			 saveFinalSeqAsFASTA=saveAsFASTA;
			 saveFinalSeqAsFASTQ=saveAsFASTQ;
		  } 
		  /*
		  barcode_saveAsFASTA=saveAsFASTA;
		  barcode_saveAsFASTQ=saveAsFASTQ;
		  primer_saveAsFASTA=saveAsFASTA;
		  primer_saveAsFASTQ=saveAsFASTQ;
		  primerCont_saveAsFASTA=saveAsFASTA;
		  primerCont_saveAsFASTQ=saveAsFASTQ;
		  bait_saveAsFASTA=saveAsFASTA;
		  bait_saveAsFASTQ=saveAsFASTQ;
		  baitBrk_saveAsFASTA=saveAsFASTA;
		  baitBrk_saveAsFASTQ=saveAsFASTQ;
		  baitArm_saveAsFASTA=saveAsFASTA;
		  baitArm_saveAsFASTQ=saveAsFASTQ; 
		  */

	}
	  
	public void setSaveHTML(boolean isSaved){ 
	    saveFinalSeqAsHTML=isSaved;
    } 
	  
	public void setSaveSeqComponent(String name, boolean saveAsFASTA, boolean saveAsFASTQ){ 
		  
		  //barcode_save=false;
		  //primer_save=false;
		  //primerCont_save=false;
		  //bait_save=false;
		  //baitBrk_save=false;
		  //baitArm_save=false;
		  
		  if(name.equalsIgnoreCase(BARCODE_NAME_DEFINITION) || name==null) {
			  //barcode_save=false;
			  //if(saveAsFASTA || saveAsFASTQ) barcode_save=true;
			  barcode_saveAsFASTA=saveAsFASTA; 
			  barcode_saveAsFASTQ=saveAsFASTQ; 	
		  }else if(name.equalsIgnoreCase(PRIMER_NAME_DEFINITION)) {		
			  //primer_save=false;
			  //if(saveAsFASTA || saveAsFASTQ) primer_save=true;
			  primer_saveAsFASTA=saveAsFASTA; 
			  primer_saveAsFASTQ=saveAsFASTQ;	
		  }else if(name.equalsIgnoreCase(PRIMERCONT_NAME_DEFINITION)) {	  		   
			  //primerCont_save=false;
			  //if(saveAsFASTA || saveAsFASTQ) primerCont_save=true;
			  primerCont_saveAsFASTA=saveAsFASTA; 
			  primerCont_saveAsFASTQ=saveAsFASTQ; 	 
		  }else if(name.equalsIgnoreCase(BAIT_NAME_DEFINITION)) {	  			
			  //bait_save=false;
			  //if(saveAsFASTA || saveAsFASTQ) bait_save=true;
			  bait_saveAsFASTA=saveAsFASTA; 
			  bait_saveAsFASTQ=saveAsFASTQ; 	    
		  }else if(name.equalsIgnoreCase(BAITBRK_NAME_DEFINITION)) {	  			
			  //baitBrk_save=false;
			  //if(saveAsFASTA || saveAsFASTQ) baitBrk_save=true;
			  baitBrk_saveAsFASTA=saveAsFASTA; 
			  baitBrk_saveAsFASTQ=saveAsFASTQ; 			  
		  }else if(name.equalsIgnoreCase(BAITARM_NAME_DEFINITION)) {
			  //baitArm_save=false;
			  //if(saveAsFASTA || saveAsFASTQ) baitArm_save=true;
			  baitArm_saveAsFASTA=saveAsFASTA; 
			  baitArm_saveAsFASTQ=saveAsFASTQ; 		 
		  }else if(name.equalsIgnoreCase("none") || name.equalsIgnoreCase("no")) {
			  //barcode_save=false;
			  //primer_save=false;
			  //primerCont_save=false;
			  //bait_save=false;
			  //baitBrk_save=false;
			  //baitArm_save=false;
			  barcode_saveAsFASTA=false;
			  barcode_saveAsFASTQ=false;
			  primer_saveAsFASTA=false;
			  primer_saveAsFASTQ=false;
			  primerCont_saveAsFASTA=false;
			  primerCont_saveAsFASTQ=false;
			  bait_saveAsFASTA=false;
			  bait_saveAsFASTQ=false;
			  baitBrk_saveAsFASTA=false;
			  baitBrk_saveAsFASTQ=false;
			  baitArm_saveAsFASTA=false;
			  baitArm_saveAsFASTQ=false;	
		  }

	}
	 
	void configRecognizer( List<SeqRocket> rockets){
	    
		SeqCompoRecognizer barcode;
		SeqCompoRecognizer primer;
		SeqCompoRecognizer primerCont;
		SeqCompoRecognizer bait;
		SeqCompoRecognizer baitBrk;
		SeqCompoRecognizer baitArm;
		
		SeqCompoRecognizer rc3together;
		/*
		SeqCompoRecognizer rc3PrimerCont;
		SeqCompoRecognizer rc3Primer;
		SeqCompoRecognizer rc3Barcode;
		SeqCompoRecognizer rc3SeqLinker;	
		*/
		  
		for(int i=0;i<rockets.size();i++){
		  
	       barcode=null;
	       primer=null;
	       primerCont=null;
		   bait=null;
	       baitBrk=null;	
	       baitArm=null;
	      
	       rc3together=null;
	      /*
	  	  rc3PrimerCont=null;
	  	  rc3Primer=null;
	  	  rc3Barcode=null;
	  	  rc3SeqLinker=null;
		  */
	       SeqCompoRecognizer recognizer;
		   String seqCompoName="";
		   rockets.get(i).saveRecognizedSeq=saveFinalSeq;
		   rockets.get(i).saveRecognizedSeqAsFASTA=saveFinalSeqAsFASTA;
		   rockets.get(i).saveRecognizedSeqAsFASTQ=saveFinalSeqAsFASTQ;
		   rockets.get(i).saveRecognizedSeqAsHTML=saveFinalSeqAsHTML;
		  
		   for(int k=0;k<rockets.get(i).seqRecognizers.size();k++){
			  recognizer=rockets.get(i).seqRecognizers.get(k);	
			  seqCompoName=rockets.get(i).seqCompoFeatures.compoNames.get(recognizer.index);
			  if(seqCompoName.equals(BARCODE_NAME_DEFINITION))
				barcode=recognizer;
			  if(seqCompoName.equals(PRIMER_NAME_DEFINITION))
				primer=recognizer;
			  if(seqCompoName.equals(PRIMERCONT_NAME_DEFINITION))
				primerCont=recognizer;
			  if(seqCompoName.equals(BAIT_NAME_DEFINITION))
			    bait=recognizer;
			  if(seqCompoName.equals(BAITBRK_NAME_DEFINITION))
			    baitBrk=recognizer;
			  if(seqCompoName.equals(BAITARM_NAME_DEFINITION))
			    baitArm=recognizer;	
			  if(seqCompoName.equals(RC3TOGETHER_NAME_DEFINITION))
				rc3together=recognizer;
			  /*
			  if(seqCompoName.equals(RC3PRIMERCONT_NAME_DEFINITION))
				rc3PrimerCont=recognizer;
			  if(seqCompoName.equals(RC3PRIMER_NAME_DEFINITION))
				rc3Primer=recognizer;
			  if(seqCompoName.equals(RC3BARCODE_NAME_DEFINITION))
				rc3Barcode=recognizer;
			  if(seqCompoName.equals(RC3SEQLINKER_NAME_DEFINITION))
				rc3SeqLinker=recognizer;
			  */
		   }
		   recognizer=null;
		  
		   if(barcode!=null){
			  barcode.minSeqLen=barcodeMinLen;
			  if(barcode.rawSeq==null && primer.rawSeq!=null){
			     barcode.seq=primer.rawSeq.substring(0,Math.min(primer.rawSeq.length(),barcodeMinLen));
			     barcode.rawSeqLength=0;
			  }else if(barcode.rawSeq!=null && barcode.rawSeq.length()<barcodeMinLen 
					&& primer.rawSeq!=null && primer.rawSeq.length()>0){			
			     barcode.seq=barcode.rawSeq+primer.rawSeq.substring(
					0, Math.min(primer.rawSeq.length(),barcodeMinLen-barcode.rawSeq.length())
			     );
			     barcode.rawSeqLength=barcode.rawSeq.length();
			  }else{
			     barcode.seq=barcode.rawSeq;
			     barcode.rawSeqLength=barcode.rawSeq.length();
			  }
			  barcode.seqLength=barcode.seq.length();
			  barcode.territoryLen=barcode.seqLength+(int) Math.ceil(barcode.seqLength*territoryLeftExtendRatio);	
			  if(barcode.rawSeqLength<=barcodeMinLen) barExactMaxStart=1;
			  barcode.exactMaxStart=barExactMaxStart;	
			
	          barcode.maxMismatchRatio=barMaxMismatchRatio;
			  barcode.maxMismatchNum=(int) Math.round(barcode.seqLength*barcode.maxMismatchRatio);
			  barcode.maxGapRatio=barMaxGapRatio;
			  barcode.maxGapNum=(int) Math.round(barcode.seqLength*barcode.maxGapRatio);	          	
	          barcode.minAlignRatio=barMinAlignRatio;
			  barcode.minAlignLen=(int) Math.ceil(barcode.seqLength*barcode.minAlignRatio);
	
			  if(barcode.minAlignLen<=12){	  
				 barcode.blastWordSize=4;
				 barcode.blastTask="blastn-short";
			  }else if(barcode.minAlignLen<=21){
				 barcode.blastWordSize=7;
				 barcode.blastTask="blastn-short";
			  }else if(barcode.minAlignLen<=50){
				 barcode.blastWordSize=11;
				 barcode.blastTask="blastn-short";
		      }else if(barcode.minAlignLen<=75){
				 barcode.blastWordSize=16;
				 barcode.blastTask="blastn";
		      }else{
			     barcode.blastWordSize=24;
				 barcode.blastTask="megablast";
			  }
			  
			  
			  //if(barcode.rawSeqLength==0){
			  if(barcode.rawSeqLength<5){
				 barcode.maxQStart=barcode.seqLength-barcode.minAlignLen+1;
				 barcode.maxSStart=2;
			  //}else if(barcode.rawSeqLength<5){
			  //	 barcode.maxQStart=1;
			  //	 barcode.maxSStart=1;
			  }else if(barcode.rawSeqLength<barcodeMinLen){
				 barcode.maxQStart=2;
				 barcode.maxSStart=2;
	          }else{			
				 barcode.maxQStart=barcode.seqLength-barcode.minAlignLen+1;
				 barcode.maxSStart=barcode.seqLength-barcode.minAlignLen+1;
			  }	   
			  barcode.leftSubForBLAST=false;
			  barcode.leftShiftForNext=false;

			  //if(primer==null) barcode.saveAsFASTAFormat=true;

			  //barcode.saveRecognizedSeq=barcode_save;
			  barcode.saveAsFASTA=barcode_saveAsFASTA;
			  barcode.saveAsFASTQ=barcode_saveAsFASTQ;
			  barcode.leftMaskSave=mask_left_barcode;	
			  barcode.leftTrimSave=trim_left_barcode;	

			  barcode.exactAlignedSeqFile=tmpDir+File.separator+barcode.seqName+"_ExactAlignedSeq";
			  barcode.seqFASTAFile=tmpDir+File.separator+"recognizer"+File.separator+barcode.seqName+".fna";	
			  SeqOperation.createFASTASeq(barcode.seq,barcode.seqName,barcode.seqFASTAFile); 

	          rockets.get(i).seqRecognizers.set(barcode.index,barcode);	
	          
			  tmpFiles.add(barcode.exactAlignedSeqFile);
			  tmpFiles.add(barcode.seqFASTAFile);

	       }
	      
		   if(primer!=null){		
			  primer.seq=primer.rawSeq;			
			  primer.seqLength=primer.seq.length();
			  if(barcode.leftShiftForNext) primer.territoryLen=primer.seqLength;
			  else primer.territoryLen= barcode.rawSeqLength+primer.seqLength;
			   
			  primer.territoryLen=primer.territoryLen+(int) Math.ceil(primer.territoryLen*territoryLeftExtendRatio);	
			  primer.exactMaxStart=primer.territoryLen-primer.seqLength+1; 			  
			  primer.maxMismatchRatio=primerMaxMismatchRatio;
			  primer.maxGapRatio=primerMaxGapRatio;		
			  primer.minAlignRatio=primerMinAlignRatio;	
			  primer.minAlignLen=(int) Math.ceil(primer.seqLength*primer.minAlignRatio);
			  primer.maxQStart=primer.seqLength-primer.minAlignLen+1;
			  primer.maxSStart=primer.territoryLen-primer.minAlignLen+1;	
			  if(primer.minAlignLen<=12){	  
				  primer.blastWordSize=4;
				  primer.blastTask="blastn-short";
			  }else if(primer.minAlignLen<=21){
				  primer.blastWordSize=7;
				  primer.blastTask="blastn-short";
			  }else if(primer.minAlignLen<=50){
				  primer.blastWordSize=11;
				  primer.blastTask="blastn-short";
			  }else if(primer.minAlignLen<=75){
				  primer.blastWordSize=16;
				  primer.blastTask="blastn";
			  }else{
				  primer.blastWordSize=24;	
				  primer.blastTask="megablast";
			  }	
			  
			  primer.leftSubForBLAST=true;
			  primer.leftShiftForNext=true;
			  if(primerCont!=null){
				  primer.trimLeftShift=primerLeftShift;
			  }else{
				  primer.trimLeftShift=0;
			  }
			  //primer.saveRecognizedSeq=primer_save;
			  primer.saveAsFASTA=primer_saveAsFASTA;
			  primer.saveAsFASTQ=primer_saveAsFASTQ;
			  primer.leftMaskSave=mask_left_primer;	
			  primer.leftTrimSave=trim_left_primer;

			  primer.seqFASTAFile=tmpDir+File.separator+"recognizer"+File.separator+primer.seqName+".fna";
			  SeqOperation.createFASTASeq(primer.seq,primer.seqName,primer.seqFASTAFile);
			  
			  rockets.get(i).seqRecognizers.set(primer.index,primer);
			  
			  tmpFiles.add(primer.seqFASTAFile);
			  
	       }	
		  
		   if(primerCont!=null){	
	          if(primerCont.rawSeq==null || primerCont.rawSeq.equals("")){
	        	  primerCont.rawSeq=bait.rawSeq.substring(0,primer.trimLeftShift);		  
	    	  }   
	          if(primer.leftShiftForNext){		  
	        	  primerCont.seq=primer.seq.substring(primer.seqLength-primer.trimLeftShift,primer.seqLength)
	        			         +primerCont.rawSeq;	
	          }else{
	        	  primerCont.seq=primer.seq+primerCont.rawSeq;
			  }
			  
	          primerCont.seqLength=primerCont.seq.length();
	          primerCont.territoryLen=primerCont.seqLength;
	          primerCont.territoryLen=primerCont.territoryLen+(int) Math.ceil(primerCont.territoryLen*territoryLeftExtendRatio);	
	          primerCont.exactMaxStart=primerCont.territoryLen-primerCont.seqLength+1;	         
	          primerCont.maxMismatchRatio=primerContMaxMismatchRatio;
	          primerCont.maxGapRatio=primerContMaxGapRatio;
	          primerCont.minAlignRatio=primerContMinAlignRatio;
	          primerCont.minAlignLen=(int) Math.ceil(primerCont.seqLength*primerCont.minAlignRatio);
			 		 
			  if(primerCont.minAlignLen<=12){	  
				  primerCont.blastWordSize=4;
				  primerCont.blastTask="blastn-short";
			  }else if(primerCont.minAlignLen<=21){
				  primerCont.blastWordSize=7;
				  primerCont.blastTask="blastn-short";
			  }else if(primerCont.minAlignLen<=50){
				  primerCont.blastWordSize=11;
				  primerCont.blastTask="blastn-short";
			  }else if(primerCont.minAlignLen<=75){
				  primerCont.blastWordSize=16;
				  primerCont.blastTask="blastn";
			  }else{
				  primerCont.blastWordSize=24;	
				  primerCont.blastTask="megablast";
			  }	
			 
			  primerCont.maxQStart=primerCont.seqLength-primerCont.minAlignLen+1;
			  primerCont.maxSStart=primerCont.territoryLen-primerCont.minAlignLen+1;	 		  
			  primerCont.leftSubForBLAST=true;
			  primerCont.leftShiftForNext=false;
			  primerCont.trimLeftShift=0;   
			  
			  //if(bait!=null){
				//  primerCont.trimLeftShift=primerCont.rawSeq.length();
	          //}else{
	        	//  primerCont.trimLeftShift=0;          	  
			  //}
	        	  
	    	  //primerCont.saveRecognizedSeq=primerCont_save;
	    	  primerCont.saveAsFASTA=primerCont_saveAsFASTA;
	    	  primerCont.saveAsFASTQ=primerCont_saveAsFASTQ;
	    	  primerCont.leftMaskSave=mask_left_primerCont;
	    	  primerCont.leftTrimSave=trim_left_primerCont;

			  primerCont.seqFASTAFile=tmpDir+File.separator+"recognizer"+File.separator+primerCont.seqName+".fna";
			  SeqOperation.createFASTASeq(primerCont.seq,primerCont.seqName,primerCont.seqFASTAFile);
			 
			  rockets.get(i).seqRecognizers.set(primerCont.index,primerCont);
			  
			  tmpFiles.add(primerCont.seqFASTAFile);
			  
	       }	

	       if(bait!=null){		
			  bait.seq=bait.rawSeq;	
			  //if(primerCont.leftShiftForNext){	 
			    //bait.seq=primerCont.seq.substring(primerCont.seqLength-primerCont.trimLeftShift,primerCont.seqLength)+bait.seq;	
			  //}else 
			  if(primer.leftShiftForNext){	 
		        bait.seq=primer.seq.substring(primer.seqLength-primer.trimLeftShift,primer.seqLength)+bait.seq;	
		      }else if(!barcode.leftShiftForNext){
		        bait.seq=primer.seq+bait.seq;
		      }else{
		        bait.seq=barcode.rawSeq+primer.seq+bait.seq;
		      }
		      bait.seqLength=bait.seq.length();
	          bait.territoryLen=bait.seqLength+(int) Math.ceil(bait.seqLength*territoryLeftExtendRatio);	
			  bait.exactMaxStart=bait.territoryLen-bait.seqLength+1;		          	  
			  bait.maxMismatchRatio=baitMaxMismatchRatio;
			  bait.maxGapRatio=baitMaxGapRatio;
			  bait.minAlignRatio=baitMinAlignRatio;	
			  bait.minAlignLen=primerCont.minAlignLen;
			  bait.maxQStart=primerCont.maxQStart;
		      bait.maxSStart=primerCont.maxSStart;	
	          bait.blastWordSize=primerCont.blastWordSize;		  
		      
		      bait.leftSubForBLAST=false;
			  if(bait.territoryLen<100)  bait.leftSubForBLAST=true;
			  
		      bait.leftShiftForNext=true;		  
			  if(baitBrk!=null || baitArm!=null ){
			    bait.trimLeftShift=0;
			  }else{
			    bait.trimLeftShift=baitLeftShift;		 
			  }
			  //bait.saveRecognizedSeq=bait_save;
			  bait.saveAsFASTA=bait_saveAsFASTA;
			  bait.saveAsFASTQ=bait_saveAsFASTQ;
			  bait.leftMaskSave=mask_left_bait;
			  bait.leftTrimSave=trim_left_bait;	

		      bait.seqFASTAFile=tmpDir+File.separator+"recognizer"+File.separator+bait.seqName+".fna";	
		      SeqOperation.createFASTASeq(bait.seq,bait.seqName,bait.seqFASTAFile);
			  
			  rockets.get(i).seqRecognizers.set(bait.index,bait);
			  
			  tmpFiles.add(bait.seqFASTAFile);
			 
	       }
		  
		   if(baitBrk!=null){		
			  baitBrk.seq=baitBrk.rawSeq;
		      baitBrk.seqLength=baitBrk.seq.length();
	          baitBrk.territoryLen=baitBrk.seqLength+(int) Math.ceil(baitBrk.seqLength*territoryLeftExtendRatio);	
			  baitBrk.exactMaxStart=baitBrk.territoryLen-baitBrk.seqLength+1;	         		  
			  baitBrk.maxMismatchRatio=baitBrkMaxMismatchRatio;
			  baitBrk.maxGapRatio=baitBrkMaxGapRatio;
			  baitBrk.minAlignRatio=baitBrkMinAlignRatio;
			  baitBrk.minAlignLen=baitBrkMinAlign;		
			  baitBrk.maxQStart=baitBrk.seqLength-baitBrk.minAlignLen+1;
		      baitBrk.maxSStart=baitBrk.territoryLen-baitBrk.minAlignLen+1;
			  if(baitBrk.minAlignLen<=12)	  
				baitBrk.blastWordSize=4;
			  else if(baitBrk.minAlignLen<=21)
				baitBrk.blastWordSize=7;
			  else if(baitBrk.minAlignLen<=33)
				baitBrk.blastWordSize=11;
			  else
			    baitBrk.blastWordSize=16;
			
			  baitBrk.leftSubForBLAST=false;
			  if(baitBrk.territoryLen<100)  baitBrk.leftSubForBLAST=true;
			  
			  baitBrk.leftShiftForNext=true;
			  if(baitArm!=null){
				 baitBrk.trimLeftShift=0;
			  }else{
				 baitBrk.trimLeftShift=baitBrkLeftShift;
	          }
			  //baitBrk.saveRecognizedSeq=baitBrk_save;
			  baitBrk.saveAsFASTA=baitBrk_saveAsFASTA;
			  baitBrk.saveAsFASTQ=baitBrk_saveAsFASTQ;
			  baitBrk.leftMaskSave=mask_left_baitBrk;
			  baitBrk.leftTrimSave=trim_left_baitBrk;

		      baitBrk.seqFASTAFile=tmpDir+File.separator+"recognizer"+File.separator+baitBrk.seqName+".fna";	
		      SeqOperation.createFASTASeq(baitBrk.seq,baitBrk.seqName,baitBrk.seqFASTAFile);
			  
			  rockets.get(i).seqRecognizers.set(baitBrk.index,baitBrk);
			  
			  tmpFiles.add(baitBrk.seqFASTAFile);			 
	       }
		  
		   if(baitArm!=null){		
			  baitArm.seq=baitArm.rawSeq;
		      baitArm.seqLength=baitArm.seq.length();
	          baitArm.territoryLen=baitArm.seqLength+(int) Math.ceil(baitArm.seqLength*territoryLeftExtendRatio);	
			  baitArm.exactMaxStart=baitArm.territoryLen-baitArm.seqLength+1;	         	  
			  baitArm.maxMismatchRatio=baitArmMaxMismatchRatio;
			  baitArm.maxGapRatio=baitArmMaxGapRatio;
			  baitArm.minAlignRatio=baitArmMinAlignRatio;	
			  baitArm.minAlignLen=baitArmMinAlign;			
			  baitArm.maxQStart=baitArm.seqLength-baitArm.minAlignLen+1;
		      baitArm.maxSStart=baitArm.territoryLen-baitArm.minAlignLen+1;	
	          if(baitArm.minAlignLen<=12)	  
				baitArm.blastWordSize=4;
			  else if(baitArm.minAlignLen<=21)
				baitArm.blastWordSize=7;
			  else if(baitArm.minAlignLen<=33)
				baitArm.blastWordSize=11;
			  else
			    baitArm.blastWordSize=16;		
				
		      baitArm.leftSubForBLAST=false;
			  if(baitArm.territoryLen<100)  baitArm.leftSubForBLAST=true;
			  
		      baitArm.leftSubForBLAST=true;
		      baitArm.leftShiftForNext=true;
			  baitArm.trimLeftShift=baitArmLeftShift;
		      //baitArm.saveRecognizedSeq=baitArm_save;
			  baitArm.saveAsFASTA=baitArm_saveAsFASTA;
			  baitArm.saveAsFASTQ=baitArm_saveAsFASTQ;
			  baitArm.leftMaskSave=mask_left_baitArm;
		      baitArm.leftTrimSave=trim_left_baitArm;
			
		      baitArm.seqFASTAFile=tmpDir+File.separator+"recognizer"+File.separator+baitArm.seqName+".fna";	
		      SeqOperation.createFASTASeq(baitArm.seq,baitArm.seqName,baitArm.seqFASTAFile);
			  
			  rockets.get(i).seqRecognizers.set(baitArm.index,baitArm);
			  
			  tmpFiles.add(baitArm.seqFASTAFile);			 
	       }
		  
		   if(rc3together!=null){	
			  rc3together.seq=rc3together.rawSeq;		  
			  rc3together.seqLength=rc3together.seq.length();
			  rc3together.territoryLen=rc3together.seqLength;
			  rc3together.territoryLen=rc3together.territoryLen+(int) Math.ceil(rc3together.territoryLen*territoryRightExtendRatio);	
			  rc3together.territoryPercent=rightTerritoryPercent;
			  rc3together.maxMismatchRatio=rc3MaxMismatchRatio;			  
			  rc3together.maxGapRatio=rc3MaxGapRatio;
			  rc3together.minAlignRatio=rc3MinAlignRatio;
			  rc3together.minAlignLen=rc3MinAlignLen;
			  rc3together.maxQStart=(int) Math.ceil(rc3together.seqLength*rc3MaxTStartRatio);	 
			  if(rc3together.minAlignLen<=12){	  
				  rc3together.blastWordSize=4;
				  rc3together.blastTask="blastn-short";
			  }else if(rc3together.minAlignLen<=21){
				  rc3together.blastWordSize=7;
				  rc3together.blastTask="blastn-short";
			  }else if(rc3together.minAlignLen<=50){
				  rc3together.blastWordSize=11;
				  rc3together.blastTask="blastn-short";
			  }else if(rc3together.minAlignLen<=75){
				  rc3together.blastWordSize=16;
				  rc3together.blastTask="blastn";
			  }else{
				  rc3together.blastWordSize=24;	
				  rc3together.blastTask="megablast";
			  }	
			 
	 		  
			  rc3together.rightSubForBLAST=false;
			  rc3together.rightShiftForNext=false;

			  rc3together.trimRightShift=0;
			  //rc3together.saveRecognizedSeq=false;
			  rc3together.saveAsFASTA=false;
			  rc3together.saveAsFASTQ=false;
			  rc3together.rightMaskSave=mask_right_rc3together;        
			  rc3together.rightTrimSave=trim_right_rc3together;  
			 		  
			  rc3together.seqFASTAFile=tmpDir+File.separator+"recognizer"+File.separator+rc3together.seqName+".fna";
			  SeqOperation.createFASTASeq(rc3together.seq,rc3together.seqName,rc3together.seqFASTAFile);
			 
			  rockets.get(i).seqRecognizers.set(rc3together.index,rc3together);
			  
			  tmpFiles.add(rc3together.seqFASTAFile);			  
	       }	
		 	  
		   barcode=null;
		   primer=null;
		   primerCont=null;
		   bait=null; 
		   baitBrk=null;
		   baitArm=null;
		  
		   rc3together=null;
		
		}	
	   
	}  	 

}
