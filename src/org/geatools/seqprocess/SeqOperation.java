package org.geatools.seqprocess;
/*******************************************************************************
 *  ========================================================================
 *  GEATools : a free Genomic Event Analysis Tools for the Java(tm) platform
 *  ========================================================================
 *
 *  (C) Copyright 2016, by Qi Wang and Contributors.
 *
 *  This file is part of GEATools.
 *
 *  GEATools is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License as published by
 *  the Free Software Foundation; either version 2.1 of the License, or
 *  (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful, but
 *  WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
 *  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>
 *
 *  [Oracle and Java are registered trademarks of Oracle and/or its affiliates. 
 *  Other names may be trademarks of their respective owners.]
 *
 *  (C) Copyright 2000-2014, by Original Author and Contributors.
 *
 *  Original Author: Qi Wang;
 *  Contributor(s): 
 *  Changes (from 1-Jan-2016)
 *  ---------------------------------------------------------------------------
 *******************************************************************************/
import java.io.*;  
import java.util.ArrayList; 
import java.util.List;

import org.geatools.data.structure.ChrSite;
import org.geatools.data.structure.SeqChrAlignSite;
import org.geatools.data.structure.SeqInfo;
import org.geatools.operation.FileOperate;

public class SeqOperation{

 public static int nameColIdx=0;
 public static int nameStartRowIdx=1;
 public static final String SEQTYPE_SINGLEEND="SingleEnd";
 public static final String SEQTYPE_PAIREND="PairEnd";
 
 static int splitStep=500000;
 static String tmpDir;
 static List<String> tmpFiles; 
 public SeqOperation(){

 }
 
 public void setTmpDir(String dir){
	tmpDir=dir;
 }
 public static int runBLAST2Seq(String blastCMD){
	   try {
	        Runtime rt = Runtime.getRuntime();        
	        Process pr = rt.exec(blastCMD);
	 
	        BufferedReader input = new BufferedReader(
	        		new InputStreamReader(pr.getInputStream()));
	        String line=null;
	        while((line=input.readLine()) != null) {
	            System.out.println(line);
	        }
	 
	        int exitVal = pr.waitFor();
			if(exitVal>0) System.out.println("Error: Exited with error code "+exitVal);
			pr=null;
			
			return exitVal;
	 
	    } catch(Exception e) {
	        System.out.println(e.toString());
	        e.printStackTrace();
	        
	        return 1;
	    }
	    
 }
 
 public static List<SeqInfo> getSeqObj(String seqFile){
     
	 List<SeqInfo> seqObj=new ArrayList<SeqInfo>();
	 if(isFASTASeq(seqFile)){		
		 seqObj=getFASTASeqObj(seqFile);		 
	 }else if(isFASTQSeq(seqFile)){		 
		 seqObj=getFASTQSeqObj(seqFile);		 
	 }
	 
	 return seqObj;
		
 }
 
 public static  List<SeqInfo> getFASTASeqObj(String seqFile){
 
    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	SeqInfo seqObj;
	int seqNum=0;
	try{    
        BufferedReader br;              
        br = new BufferedReader(new FileReader(seqFile));
	    String line;
		String seqIdentifier;
		String seqName;
		String seqLine;
		String [] itemSplited;
			
		seqNum=0;
		line = br.readLine();				
		while(true){           
	       if (line == null) break;
	       line=line.trim();
		   if(line.indexOf(">")==0){
		     
			 seqIdentifier=line.substring(1,line.length());			 
			 itemSplited=seqIdentifier.split("\\s+");
			 seqName=itemSplited[0].trim();	
			 
			 line=br.readLine();
			 seqLine="";
			 if (line == null) break;
			 while(line.indexOf(">")<0){
		      seqLine = seqLine+line.trim();
			  line=br.readLine();
			  if (line == null) break;
			 }
			 seqLine=seqLine.replaceAll("N","n");
             seqObj=new SeqInfo();
             seqObj.seqIdentifier=seqIdentifier;
             seqObj.seqName=seqName;
			 seqObj.seqLength=seqLine.length();
             seqObj.seq=seqLine.toUpperCase();
             seqObjList.add(seqObj);
             seqObj=null;
             
             seqNum=seqNum+1;
             
		   }
			 
		}           		   
		br.close();
	
	}catch(IOException e){
        System.out.println(e);
    }
	
	return seqObjList;
   
 }
 
 public static List<SeqInfo> getFASTASeqObj(String seqFile, int start, int end){
	 
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq=new SeqInfo();
		int seqNum=0;
		
		try{    
	       	BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier="";
			String seqLine = "";
			String seqName="";			
			String [] itemSplited;
			seqNum=0;
			line = br.readLine();	
			while(true){  
	         		
		       if (line == null) break;
		       if(seqNum>=end) break;
		       line=line.trim();
			   if(line.indexOf(">")==0){			    
				 seqIdentifier=line.substring(1,line.length());				 
				 line=br.readLine();				
				 if (line == null) break;
				 seqLine="";
				 while(line.indexOf(">")<0){
			       seqLine = seqLine+line.trim();
				   line=br.readLine();
				   if (line == null) break;
				 }
				 
				 if(seqNum>=start && seqNum<end){					   
					 itemSplited=seqIdentifier.split("\\s+");
					 seqName=itemSplited[0].trim();					
					 seqLine=seqLine.replaceAll("N","n");
					 perSeq=new SeqInfo();				  
					 perSeq.seqIdentifier=seqIdentifier;
					 perSeq.seqName=seqName;
					 perSeq.seqLength=seqLine.length();
					 perSeq.seq=seqLine.toUpperCase();
			         seqObjList.add(perSeq);
			         perSeq=null;	
				 }
				 
				 seqNum=seqNum+1;
	             
			   }
			}           		   
			br.close();
		
	    	line=null;
	    	seqIdentifier=null;		
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
 }
 
 public static List<SeqInfo> getFASTQSeqObj(String seqFile){
	 
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq=new SeqInfo();
		int seqNum=0;
		
		try{    
	       	BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;	
			String seqLine;
			String seqName;		
			String seqQualityLine;
			String [] itemSplited;
				
			while(true){  			   
			   
	           line = br.readLine();			
		       if (line == null) break;
			   line=line.trim();
	           if(line.indexOf("@")!=0) {
			      System.out.println("Error in reading fastq");
				  break;
			   }			   
	      	   seqIdentifier=line.substring(1,line.length());				 
			   itemSplited=seqIdentifier.split("\\s+");
			   seqName=itemSplited[0].trim();
			     
			   line=br.readLine();			 
			   if (line == null) break;
			   seqLine=line.trim();			 
			   //seqLine=seqLine.replaceAll("N","n");
			   
			   line=br.readLine();			  
			   if (line == null) break;
			   line=line.trim();
			   if(line.indexOf("+")!=0) {
			      System.out.println("Error in reading fastq");
				  break;
			   }		      
	           
	           line=br.readLine();			 
			   if (line == null) break;
			   seqQualityLine=line.trim();
			   
			   perSeq=new SeqInfo();
			   perSeq.seqIdentifier=seqIdentifier;
			   perSeq.seqName=seqName;
			   perSeq.seqLength=seqLine.length();
			   perSeq.seq=seqLine;
			   perSeq.seqQualityEncode=seqQualityLine;
			   seqObjList.add(perSeq);
			   perSeq=null;
			   
			   seqNum=seqNum+1;
		
			}           		   
			br.close();
		
	    	line=null;
	    	seqIdentifier=null;		
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
 }
 
 public static List<SeqInfo> getFASTQSeqObj(String seqFile, int start, int end){
	 
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq=new SeqInfo();
		int seqNum=0;
		
		try{    
	       	BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;		
			String seqLine;
			String seqName;		
			String seqQualityLine;
			String [] itemSplited;
				
			while(true){  
			   
		       if(seqNum>=end) break;
		       
	           line = br.readLine();			
		       if (line == null) break;
		       line=line.trim();
		       if(line.indexOf("@")!=0) {
			      System.out.println("Error in reading fastq");
				  break;
			   }			   
		       seqIdentifier=line.substring(1,line.length());
			     
			   line=br.readLine();			 
			   if (line == null) break;
			   seqLine=line.trim();			 			   
			   
			   line=br.readLine();			  
			   if (line == null) break;
			   if(line.indexOf("+")!=0) {
			      System.out.println("Error in reading fastq");
				  break;
			   }	
	                 
	           line=br.readLine();			 
			   if (line == null) break;
			   seqQualityLine=line;
			   
			   if(seqNum>=start && seqNum<end){
				   				 
				   itemSplited=seqIdentifier.split("\\s+");
				   seqName=itemSplited[0].trim();				
				   seqLine=seqLine.replaceAll("N","n");
				   
				   perSeq=new SeqInfo();
				   perSeq.seqIdentifier=seqIdentifier;
				   perSeq.seqName=seqName;
				   perSeq.seqLength=seqLine.length();
				   perSeq.seq=seqLine.trim();
				   perSeq.seqQualityEncode=seqQualityLine.trim();
				   seqObjList.add(perSeq);
				   perSeq=null;
			   }
			   
			   seqNum=seqNum+1;
		
			}           		   
			br.close();
		
	    	line=null;
	    	seqIdentifier=null;	
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
 }
 
 public static String convertFASTQ2FASTA(String inSeqFile){
	  
	    List<SeqInfo> seqObjList;	
	    tmpFiles=new ArrayList<String>();
	    int total=SeqOperation.getFASTQSeqNum(inSeqFile);
	    String outSeqFile=inSeqFile.substring(0,inSeqFile.lastIndexOf("."))+".fasta";	
		try{
		  BufferedWriter writer=null;		   
		  writer=new BufferedWriter(new FileWriter(outSeqFile));	
		  int stopPos=0;
		  for(int f=0;f<total;f=f+splitStep){	
	         String seqIdentifier;
	         String seq;	
	         stopPos=Math.min(total, (f+splitStep));
	         seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;					
				writer.write(">"+seqIdentifier);
				writer.newLine();		
				writer.write(seq);
				writer.newLine();

				writer.flush();			
		     }
		  }
		  writer.close();
		}catch(IOException e){
	        System.out.println(e);
	    }
		seqObjList=null;
		
		return outSeqFile;

 } 

 public static boolean isFASTASeq(List<String> seqFileList){
	 boolean isOK=false;
	 for(String seqFile:seqFileList){
		 if(isFASTASeq(seqFile)){
			 isOK=true;
			 break;
		 }
	 }
	 
	 return isOK;
 }
 public static boolean isFASTASeq(String seqFile){
	 if(seqFile==null) return false;
	 
	 File f=new File(seqFile);
	 if(!f.exists()){
		f=null; 
		return false;
	 }
	   
	 boolean isOk=false;
	 
	 String inSeqFileFormat=FileOperate.getFileFormat(seqFile);
	
     if(inSeqFileFormat.equalsIgnoreCase("fasta") 
			|| inSeqFileFormat.equalsIgnoreCase("fna") 
			|| inSeqFileFormat.equalsIgnoreCase("fa")){
		
		
    	 isOk=true;
	
	 }
	 
	 return isOk;
 }
 
 public static boolean isFASTQSeq(List<String> seqFileList){
	 boolean isOK=false;
	 for(String seqFile:seqFileList){
		 if(isFASTQSeq(seqFile)){
			 isOK=true;
			 break;
		 }
	 }
	 
	 return isOK;
 }
 
 public static boolean isFASTQSeq(String seqFile){
	 if(seqFile==null) return false;
	 
	 File f=new File(seqFile);
	 if(!f.exists()){
		f=null; 
		return false;
	 }
	 
	 boolean isOk=false;
	 
	 String inSeqFileFormat=FileOperate.getFileFormat(seqFile);
	
     if(inSeqFileFormat.equalsIgnoreCase("fastq") 
			|| inSeqFileFormat.equalsIgnoreCase("fnq") 
			|| inSeqFileFormat.equalsIgnoreCase("fq")){
		
		
    	 isOk=true;
	
	 }
	 
	 return isOk;
 }
 
 public static boolean isFASTAFormat(String seqFile){
	 if(seqFile==null) return false;	 
	 
	 boolean isOk=false;
	 
	 String inSeqFileFormat=FileOperate.getFileFormat(seqFile);
	
     if(inSeqFileFormat.equalsIgnoreCase("fasta") 
			|| inSeqFileFormat.equalsIgnoreCase("fna") 
			|| inSeqFileFormat.equalsIgnoreCase("fa")){
		
		
    	 isOk=true;
	
	 }
	 
	 return isOk;
 }
 
 public static boolean isFASTQFormat(String seqFile){
	 if(seqFile==null) return false;	 
	 
	 boolean isOk=false;
	 
	 String inSeqFileFormat=FileOperate.getFileFormat(seqFile);
	
     if(inSeqFileFormat.equalsIgnoreCase("fastq") 
			|| inSeqFileFormat.equalsIgnoreCase("fnq") 
			|| inSeqFileFormat.equalsIgnoreCase("fq")){
		
		
    	 isOk=true;
	
	 }
	 
	 return isOk;
 }
 
 public static List<String> splitSeqFile(String inSeqFile, int step, String outDir){
	     
	     List<String> splitedSeqFiles=new ArrayList<String>();		
	     
	     if(isFASTASeq(inSeqFile)){			
			splitedSeqFiles=SeqOperation.splitFASTAFile(inSeqFile, step, outDir);		
		 }else if(isFASTQSeq(inSeqFile)){			 
			splitedSeqFiles=SeqOperation.splitFASTQFile(inSeqFile, step, outDir);
		 }
	     
	     return splitedSeqFiles;
 }
 
 public static List<String> splitFASTAFile(String inSeqFile, int step, String outDir){
	    	   
		List<String> splitedFiles=new ArrayList<String>();
	    List<SeqInfo> seqObjList;	
		if(step<0) step=500000;
		if(outDir==null) 
		  outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/"
		       +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
		       +"_Splited";
		FileOperate.newFolder(outDir);
		int total=SeqOperation.getFASTASeqNum(inSeqFile);
		System.out.println("Total Seq Counts: "+total);    	
		try{
		  BufferedWriter writer=null;
		  int stopPos=0;		 
		  for(int f=0;f<total;f=f+step){	
		     
		     String outSeqFile=outDir+"/split."+step+"_"+(f+1)+"-"+(f+step)+".fna";	  
		     writer=new BufferedWriter(new FileWriter(outSeqFile));	
	         String seqIdentifier;
	         String seq;	      
	         stopPos=Math.min(total, (f+step));
	         seqObjList=getFASTASeqObj(inSeqFile,f,stopPos);	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
	 			seq=seqObjList.get(i).seq;
	 			writer.write(">"+seqIdentifier);
	 			writer.newLine();
	 			//writer.flush(); 
	 			writer.write(seq);
	 			writer.newLine();
	 			writer.flush();				
		     }
		     writer.close();
		     splitedFiles.add(outSeqFile);
		     
		     System.out.println("Split "+(f+1)+"-"+(f+step)+" OK!");	
		  }
		}catch(IOException e){
	        System.out.println(e);
	    }
		seqObjList=null;
		  
		return  splitedFiles;
 }
 
 public static List<String> splitFASTQFile(String inSeqFile, int step, String outDir){
	    	   
		List<String> splitedFiles=new ArrayList<String>();
	    List<SeqInfo> seqObjList;	
		if(step<0) step=500000;
		if(outDir==null) 
			outDir=inSeqFile.substring(0,inSeqFile.lastIndexOf("/"))+"/"
		           +inSeqFile.substring(inSeqFile.lastIndexOf("/")+1,inSeqFile.lastIndexOf("."))
		           +"_Splited";
		FileOperate.newFolder(outDir);
		int total=SeqOperation.getFASTQSeqNum(inSeqFile);
		System.out.println("Total Seq Counts: "+total);    	
		try{
		  BufferedWriter writer=null;
		  int stopPos=0;
		 
		  for(int f=0;f<total;f=f+step){	
		     
		     String outSeqFile=outDir+"/split."+step+"_"+(f+1)+"-"+(f+step)+".fastq";	  
		     writer=new BufferedWriter(new FileWriter(outSeqFile));	
	         String seqIdentifier;
	         String seq;
	         String seqQuality;
	         stopPos=Math.min(total, (f+step));
	         seqObjList=getFASTQSeqObj(inSeqFile,f,stopPos);	
	         for(int i=0;i<seqObjList.size();i++){	
	        	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;
				seqQuality=seqObjList.get(i).seqQualityEncode;
				
				writer.write("@"+seqIdentifier);
				writer.newLine();		
				writer.write(seq);
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(seqQuality);
				writer.newLine();
				writer.flush();			
		     }
		     writer.close();
		     splitedFiles.add(outSeqFile);
		     
		     System.out.println("Split "+(f+1)+"-"+(f+step)+" OK!");	
		  }
		}catch(IOException e){
	        System.out.println(e);
	    }
		seqObjList=null;
		  
		return  splitedFiles;
 }

 public static String combineSeqFile(List<String> splitedFiles, String outFile){
	 
	 if(splitedFiles==null || splitedFiles.size()==0) return null;
	 String seqFile=splitedFiles.get(0);
	 
	 if(isFASTASeq(seqFile)){		
		outFile=combineFASTAFile(splitedFiles, outFile);		 
	 }else if(isFASTQSeq(seqFile)){		 
		outFile=combineFASTQFile(splitedFiles, outFile);		 
	 }
	 
	 return outFile;
		
 }
 
 public static String combineFASTAFile(List<String> splitedFiles, String outFile){
	    
	    String seqIdentifier;
	    String seq;	
	    String fastaFile=splitedFiles.get(0);
	   
		if(outFile==null) {
			String outDir=fastaFile.substring(0,fastaFile.lastIndexOf("/"))+"/combined";
			outFile=outDir+"/"
			        +fastaFile.substring(fastaFile.lastIndexOf("/"),fastaFile.lastIndexOf("."))
			        +"_combined.fna";
			FileOperate.newFolder(outDir);
		}
		
		System.out.println("Combining.................. ");	
		try{
		  BufferedWriter writer=null;
		  writer=new BufferedWriter(new FileWriter(outFile));
		  int num=0;
		  for(String file:splitedFiles ){	     
			 List<SeqInfo> seqObjList=getSeqObj(file);	 
			 num=num+seqObjList.size();
		     for(int i=0;i<seqObjList.size();i++){		
		    	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;
				writer.write(">"+seqIdentifier);
				writer.newLine();			
				writer.write(seq);
				writer.newLine();
				writer.flush();			
		     }	     
		     seqObjList=null;
		  }
		  writer.close();
		  System.out.println("Totally combined sequences number: "+num);

		}catch(IOException e){
	        System.out.println(e);
	    }
			  
		return  outFile;
 }
 
 public static String combineFASTQFile(List<String> splitedFiles, String outFile){
		    
	    String seqIdentifier;
	    String seq;	
	    String seqQuality;
	    String fastqFile=splitedFiles.get(0);
	   
		if(outFile==null) {
			String outDir=fastqFile.substring(0,fastqFile.lastIndexOf("/"))+"/combined";
			outFile=outDir+"/"
			        +fastqFile.substring(fastqFile.lastIndexOf("/"),fastqFile.lastIndexOf("."))
			        +"_combined.fastq";
			FileOperate.newFolder(outDir);
		}
		
		System.out.println("Combining.................. ");	
		try{
		  BufferedWriter writer=null;
		  writer=new BufferedWriter(new FileWriter(outFile));
		  int num=0;
		  for(String file:splitedFiles ){	     
			 List<SeqInfo> seqObjList=getSeqObj(file);	 
			 num=num+seqObjList.size();
		     for(int i=0;i<seqObjList.size();i++){		
		    	seqIdentifier=seqObjList.get(i).seqIdentifier;
				seq=seqObjList.get(i).seq;
				seqQuality=seqObjList.get(i).seqQualityEncode;
				
				writer.write("@"+seqIdentifier);
				writer.newLine();		
				writer.write(seq);
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(seqQuality);
				writer.newLine();
				writer.flush();		
				writer.flush();			
		     }	     
		     seqObjList=null;
		  }
		  writer.close();
		  System.out.println("Totally combined sequences number: "+num);

		}catch(IOException e){
	        System.out.println(e);
	    }
			  
		return  outFile;
 }
	 
 public List<SeqInfo> checkSeq(String inSeqFile){
	    
	    List<SeqInfo> seqObjList = null;

		if(isFASTASeq(inSeqFile)){		  
			seqObjList=checkFASTASeq(inSeqFile);			
		}else if(isFASTQSeq(inSeqFile)){		  
			seqObjList=checkFASTQSeq(inSeqFile);			
		}
		
		return seqObjList;
 }
 
 static List<SeqInfo> checkFASTASeq(String seqFile, String outSeqFile){
	  
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq=new SeqInfo();
		int seqNum=0;
		//ArrayList<String> uniSeqName=new ArrayList<String> ();
		try{    
	        BufferedWriter writer=null;
	        writer=new BufferedWriter(new FileWriter(outSeqFile));
			BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;
			String seqLine;
			String seqName;
			String [] itemSplited;
			//uniSeqName=new ArrayList<String> ();
			seqObjList=new ArrayList<SeqInfo>();
			seqNum=0;
			line = br.readLine();				
			while(true){           
		       if (line == null) break;
			   if(line.trim().indexOf(">")==0){
			     seqNum=seqNum+1;
			     seqIdentifier=line.substring(1,line.length());				 
				 itemSplited=seqIdentifier.split("\\s+");
				 seqName=itemSplited[0].trim();
				  
				 line=br.readLine();
				 seqLine="";
				 if (line == null) break;
				 while(line.indexOf(">")<0){
			      seqLine = seqLine+line.trim();
				  line=br.readLine();
				  if (line == null) break;
				 }
				 seqLine=seqLine.replaceAll("N","n");

				//uniSeqName.add(seqName);
				 writer.write(">"+seqIdentifier);
				 writer.newLine();
				 writer.write(seqLine);
				 writer.newLine();
				 writer.flush();  
				  
				 perSeq=new SeqInfo();
				 perSeq.seqIdentifier=seqIdentifier;
				 perSeq.seqName=seqName;
				 perSeq.seqLength=seqLine.length();
				 perSeq.seq=seqLine;
				 seqObjList.add(perSeq);
				 perSeq=null;
	             
			   }
				 
			}           		   
			br.close();
			writer.close();
			
			line=null;
			seqIdentifier=null;	
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
 }
	 
 static List<SeqInfo> checkFASTQSeq(String seqFile, String outSeqFile){
	  
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq=new SeqInfo();
		try{    
	        BufferedWriter writer=null;
	        writer=new BufferedWriter(new FileWriter(outSeqFile));
			BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;		
			String seqLine;
			String seqName;		
			String seqQualityLine;
			String [] itemSplited;
		
			seqObjList=new ArrayList<SeqInfo>();
			
			while(true){    
			   line = br.readLine();	// seq name		
		       if (line == null) break;
			   line=line.trim();
	           if(line.indexOf("@")!=0) {
			      System.out.println("Error in reading fastq");
				  break;
			   }			   
	           seqIdentifier=line.substring(1,line.length());			 
			   itemSplited=seqIdentifier.split("\\s+");
			   seqName=itemSplited[0].trim();
			     
			   line=br.readLine();		// seq line	 
			   if (line == null) break;
			   seqLine=line.trim();			 
			   seqLine=seqLine.replaceAll("N","n");
			   
			   line=br.readLine();		// '+' tag	  
			   if (line == null) break;
			   line=line.trim();
			   if(line.indexOf("+")!=0) {
			      System.out.println("Error in reading fastq");
				  break;
			   }	   
	           
	           line=br.readLine();		// seq quality	 
			   if (line == null) break;
			   seqQualityLine=line.trim();
				
			   writer.write("@"+seqIdentifier);
			   writer.newLine();
			   writer.write(seqLine);
			   writer.newLine();
			   writer.write("+");
			   writer.newLine();
			   writer.write(seqQualityLine);
			   writer.newLine();
			   writer.flush();  
				  
			   perSeq=new SeqInfo();
			   perSeq.seqIdentifier=seqIdentifier;
			   perSeq.seqName=seqName;
			   perSeq.seqLength=seqLine.length();
			   perSeq.seq=seqLine;
			   perSeq.seqQualityEncode=seqQualityLine;
			   seqObjList.add(perSeq);
			   perSeq=null;
				 
			}           		   
			br.close();
			writer.close();
			
			line=null;
			seqIdentifier=null;
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
  }
	 
  static List<SeqInfo> checkFASTASeq(String seqFile){
	  
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq=new SeqInfo();
		int seqNum=0;

		try{
			
			BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;
			String seqLine;
			String seqName;
			String [] itemSplited;
			seqObjList=new ArrayList<SeqInfo>();
			seqNum=0;
			line = br.readLine();				
			while(true){           
		       if (line == null) break;
			   if(line.indexOf(">")==0){
			     seqNum=seqNum+1;
			     seqIdentifier=line.substring(1,line.length());		 
				 itemSplited=seqIdentifier.split("\\s+");
				 seqName=itemSplited[0].trim();
				 
				 line=br.readLine();
				 seqLine="";
				 if (line == null) break;
				 while(line.indexOf(">")<0){
			      seqLine = seqLine+line.trim();
				  line=br.readLine();
				  if (line == null) break;
				 }
				 seqLine=seqLine.replaceAll("N","n");
			  
				 perSeq=new SeqInfo();
				 perSeq.seqIdentifier=seqIdentifier;
				 perSeq.seqName=seqName;
				 perSeq.seqLength=seqLine.length();
				 perSeq.seq=seqLine;
				 seqObjList.add(perSeq);
				 perSeq=null;
	             
			   }
				 
			}           		   
			br.close();
		
			line=null;
			seqIdentifier=null;	
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
  }
	 
  static List<SeqInfo> checkFASTQSeq(String seqFile){
	  
	    List<SeqInfo> seqObjList=new ArrayList<SeqInfo>();
	    SeqInfo perSeq=new SeqInfo();
			
		try{    
	       	BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;
			String seqIdentifier;		
			String seqLine;
			String seqName;		
			String seqQualityLine;
			String [] itemSplited;
		
			seqObjList=new ArrayList<SeqInfo>();
			
			while(true){  
			  
	           line = br.readLine();		// seq name	
		       if (line == null) break;
			   line=line.trim();
	           if(line.indexOf("@")!=0) {
			      System.out.println("Error in reading fastq");
				  break;
			   }			   
	           seqIdentifier=line.substring(1,line.length());				 
			   itemSplited=seqIdentifier.split("\\s+");
			   seqName=itemSplited[0].trim();
			     
			   line=br.readLine();		// seq line	 
			   if (line == null) break;
			   seqLine=line.trim();			 
			   seqLine=seqLine.replaceAll("N","n");
			   
			   line=br.readLine();		// '+' tag	  
			   if (line == null) break;
			   line=line.trim();
			   if(line.indexOf("+")!=0) {
			      System.out.println("Error in reading fastq");
				  break;
			   }	
	          	           
	           line=br.readLine();		// seq quality	 
			   if (line == null) break;
			   seqQualityLine=line.trim();
			   
			   perSeq=new SeqInfo();
			   perSeq.seqIdentifier=seqIdentifier;
			   perSeq.seqName=seqName;
			   perSeq.seqLength=seqLine.length();
			   perSeq.seq=seqLine;
			   perSeq.seqQualityEncode=seqQualityLine;
			   seqObjList.add(perSeq);
			   perSeq=null;
		
			}           		   
			br.close();
		
	    	line=null;
	    	seqIdentifier=null;		
			seqLine=null;
			seqName=null;
			itemSplited=null;
		}catch(IOException e){
	        System.out.println(e);
	    }

	    return seqObjList; 
 }
  
 static boolean checkPairEndSeq(List<SeqInfo> seqObjList, List<SeqInfo> seqObjList2){
		 boolean isOK=true;
		 if(seqObjList.size()!=seqObjList2.size()) return false;
		 
		 for(int i=0;i<seqObjList.size();i++){
			if(seqObjList.get(i).seqName.trim()!=seqObjList2.get(i).seqName.trim()){ 
				isOK=false;
				break;
			}
		 }
		 return isOK;
 }
 
 public static int getSeqNum(List<SeqInfo> seqObjList){
  
    return seqObjList.size();	
	
 }
 
 public static int getSeqNum(String seqFile){
	 
	 int seqNum=0;

	 if(isFASTASeq(seqFile)){		
		 seqNum=getFASTASeqNum(seqFile);		 
	 }else if(isFASTQSeq(seqFile)){		 
		 seqNum=getFASTQSeqNum(seqFile);		 
	 }
	 
	 return seqNum;
		
 }
 
 public static int getFASTASeqNum(String seqFile){
  
    int seqNum=0;
	try{    
		BufferedReader br;              
        br = new BufferedReader(new FileReader(seqFile));
	    String line;		
		seqNum=0;						
		while(true){
           line = br.readLine();
	       if (line == null) break;
		   if(line.trim().indexOf(">")==0){
		     seqNum=seqNum+1;			
		   }			 
		}           		   
		br.close();		
	}catch(IOException e){
        System.out.println(e);
    }   
	
    return seqNum;	  
 }
 
 public static int getFASTQSeqNum(String seqFile){
	  
	    int seqNum=0;
		try{    
			BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;		
			seqNum=0;						
			while(true){
	           line = br.readLine(); // seq name line
		       if (line == null) break;
			   if(line.trim().indexOf("@")==0){
			     seqNum=seqNum+1;
			     line = br.readLine(); // seq line
				 line = br.readLine(); // tag '+' line
				 line = br.readLine(); // seq quality line		
			   }			 
			}           		   
			br.close();		
		}catch(IOException e){
	        System.out.println(e);
	    }   
		
	    return seqNum;	  
 }
 
 public static boolean extratSeqName(String seqFile, String outFile){
	    
	    boolean isOK=true;
	
	   	try{
		  if(outFile==null) outFile=seqFile+".seqName";
		  if(isFASTASeq(seqFile)){			  
			 isOK=getFASTASeqName(seqFile,outFile);			 
		  }else if(isFASTQSeq(seqFile)){			  
			 isOK=getFASTQSeqName(seqFile,outFile);			 
		  }		  
			
		}catch(Exception e){
		  isOK=false;
		  System.out.println(e);
		}
		
	   	return isOK;
	
 }
 
 public static boolean getFASTASeqName(String seqFile, String outFile){
        boolean isOK=true;
		try{    
			BufferedReader br;              
	        br = new BufferedReader(new FileReader(seqFile));
		    String line;		
		    BufferedWriter writer=null;
		    writer=new BufferedWriter(new FileWriter(outFile));
			String seqName;
			String [] itemSplited;
			while(true){
	           line = br.readLine();
		       if (line == null) break;
			   if(line.trim().indexOf(">")==0){
				 itemSplited=line.split("\\s+");
				 seqName=itemSplited[0].trim();
				 seqName=seqName.substring(line.indexOf(">")+1,seqName.length());
				 writer.write(seqName);
				 writer.newLine();
				 writer.flush();
			   }			 
			}           		   
			br.close();
			writer.close();
		}catch(IOException e){
			isOK=false;
	        System.out.println(e);
	    }   
		
	    return isOK;
 }
	 
 public static boolean getFASTQSeqName(String seqFile,String outFile){
	        boolean isOK=true;
			try{    
				BufferedReader br;              
		        br = new BufferedReader(new FileReader(seqFile));
			    String line;		
			    BufferedWriter writer=null;
			    writer=new BufferedWriter(new FileWriter(outFile));
				String seqName;
				String [] itemSplited;
			
			    while(true){		         
			       line = br.readLine(); // seq name line
			       if (line == null) break;
			       line=line.trim();
				   if(line.indexOf("@")==0){
					 itemSplited=line.split("\\s+");
					 seqName=itemSplited[0].trim();
					 seqName=seqName.substring(seqName.indexOf("@")+1,seqName.length());
					 writer.write(seqName);
					 writer.newLine();
					 writer.flush();
					 line = br.readLine(); // seq line
					 line = br.readLine(); // tag '+' line
					 line = br.readLine(); // seq quality line					
				   }			 
				}           		   
				br.close();
				writer.close();
			}catch(IOException e){
				isOK=false;
		        System.out.println(e);
		    } 
			
			return isOK;
		
  }
 
  public static List<String> getRowName(String nameFile){
  
    List<String> seqName=new ArrayList<String>();
	List<ArrayList <String>> nameList=FileOperate.getMatrixFromFile(nameFile);	
	String name="";
	for(int i=nameStartRowIdx; i<nameList.size();i++){
	  name=nameList.get(i).get(nameColIdx);
	  if(!seqName.contains(name)){
	    seqName.add(name);	   
	  }else{	    
	    //System.out.println("Repeated Seq: "+nameList.get(i));
	  }
	}
	nameList=null;
	System.out.println("Got Unique names: "+seqName.size() +" From ["+nameFile+"]");
	
    return seqName;
 }
 
 public static List<String> getRowName(String nameFile, int nameColIdx){
	  
	    List<String> seqName=new ArrayList<String>();
		List<ArrayList <String>> nameList=FileOperate.getMatrixFromFile(nameFile);	
		String name="";
		for(int i=nameStartRowIdx; i<nameList.size();i++){
		  name=nameList.get(i).get(nameColIdx);
		  if(!seqName.contains(name)){
		    seqName.add(name);	   
		  }else{	    
		    //System.out.println("Repeated name: "+nameList.get(i));
		  }
		}
		
		nameList=null;
		System.out.println("Got Unique names: "+seqName.size() +" From ["+nameFile+"]");
		
	    return seqName;
 }
 
 
 static List<String> getRowName(List <String> nameList,int start, int end){
	  
	    List<String> names=new ArrayList<String>();			
		String name="";
		end=Math.min(end, nameList.size());
		for(int i=start; i<end;i++){
		  name=nameList.get(i);
		  if(!names.contains(name)){
			  names.add(name);	   
		  }else{	    
		    System.out.println("Repeated name: "+nameList.get(i));
		  }
		}
		
	    return names;
 }
 
 
 public static List<SeqInfo> excludeSeq(String seqFile,List<String> excludedSeqNameList){
	    
	    List<SeqInfo> seqObjList = null;

	   	try{
			
	   	  if(isFASTASeq(seqFile)){
	   		 seqObjList=getFASTASeqObj(seqFile);
		  }else if(isFASTQSeq(seqFile)){
			 seqObjList=getFASTQSeqObj(seqFile);
		  }
		  List<String> seqNameList=new ArrayList<String>();
		  for(SeqInfo seqObj:seqObjList){
			seqNameList.add(seqObj.seqName);
		  }		
	
		  int i=0;
		  for(String name: excludedSeqNameList){
			 i=seqNameList.indexOf(name);
			 if(i>=0){
			  seqObjList.remove(i);
			  seqNameList.remove(i);
			 }
		  }
	
		  seqNameList=null;
			
		}catch(Exception e){
		  System.out.println(e);
		}
		
		return seqObjList;
 }
 
 public static List<SeqInfo> excludeSeq(String seqFile,String excludedSeqNameFile,
		 int excludedSeqNameColIdx){
	    
	    List<SeqInfo> seqObjList = null;

	   	try{
			
	   	  if(isFASTASeq(seqFile)){
			 seqObjList=getFASTASeqObj(seqFile);
		  }else if(isFASTQSeq(seqFile)){
			 seqObjList=getFASTQSeqObj(seqFile);
		  }
		  List<String> seqNameList=new ArrayList<String>();
		  for(SeqInfo seqObj:seqObjList){
			 seqNameList.add(seqObj.seqName);
		  }		
		  List<String> excludedSeqNameList=getRowName(excludedSeqNameFile,excludedSeqNameColIdx);
		  int i=0;
		  for(String name: excludedSeqNameList){
			 i=seqNameList.indexOf(name);
			 if(i>=0){
				seqObjList.remove(i);
				seqNameList.remove(i);
			 }
		  }
		
		  seqNameList=null;
		  excludedSeqNameList=null;
			
		}catch(Exception e){
		  System.out.println(e);
		}
		
		return seqObjList;
 }
 
 public static List<SeqInfo> excludeSeq(List<SeqInfo> seqObjList,
		 List<String> excludedSeqNameList){

	   	try{
	
		  List<String> seqNameList=new ArrayList<String>();
		  for(SeqInfo seqObj:seqObjList){
			seqNameList.add(seqObj.seqName);
		  }		
	
		  int i=0;
		  for(String name: excludedSeqNameList){
			 i=seqNameList.indexOf(name);
			 if(i>=0){
			  seqObjList.remove(i);
			  seqNameList.remove(i);
			 }
		  }
	
		  seqNameList=null;
			
		}catch(Exception e){
		  System.out.println(e);
		}
		
		return seqObjList;
 }
 
 public static List<SeqInfo> extractSubSeq(String seqFile,String subSeqNameFile,
		 int seqNameColIdx){
	    
	    List<String> totalSubSeqNames=getRowName(subSeqNameFile,seqNameColIdx);
	    int subSeqNum=totalSubSeqNames.size();
	 
	    List<SeqInfo> subSeqObjList = new ArrayList<SeqInfo>();
	    SeqInfo seqInfo;
	    int total;
	    if(isFASTASeq(seqFile))
	      total=SeqOperation.getFASTASeqNum(seqFile);
	    else if(isFASTQSeq(seqFile))
		  total=SeqOperation.getFASTQSeqNum(seqFile);
		else
		  return subSeqObjList;
	    
	    int subSeqStep=1000;
	    String seqName;
	    String seq;
		int stopPos=0;	    
	    int start;
	    int end;
	    List<String> subSeqNameList;
	   	try{
	  	    for(int f=0;f<total;f=f+splitStep){
		  		stopPos=Math.min(total, (f+splitStep));
		   		List<SeqInfo> seqObjList = null;
				if(isFASTASeq(seqFile)){
				   seqObjList=getFASTASeqObj(seqFile,f,stopPos);
				}else if(isFASTQSeq(seqFile)){
				   seqObjList=getFASTQSeqObj(seqFile,f,stopPos);
				}
				List<String> seqNameList=new ArrayList<String>();
				for(SeqInfo seqObj:seqObjList){
				  seqNameList.add(seqObj.seqName);
				}	
				
				for(start=0;start<subSeqNum; start=start+subSeqStep){			
				  end=Math.min(subSeqNum, (start+subSeqStep));
				  subSeqNameList=getRowName(totalSubSeqNames,start,end);
				  int i=0;
				  for(String name: subSeqNameList){
				    i=seqNameList.indexOf(name);
					if(i>=0){
						seqInfo=new SeqInfo();
						seqName=seqObjList.get(i).seqName;
						seq=seqObjList.get(i).seq;
						seqInfo.seqIdentifier=seqObjList.get(i).seqIdentifier;
						seqInfo.seqName=seqName;
						seqInfo.seq=seq;
						seqInfo.seqQualityEncode=seqObjList.get(i).seqQualityEncode;
						subSeqObjList.add(seqInfo);
						seqInfo=null;
					}
				  }
				  System.out.println((f+1)+"."+(start+1)+"-"+end+" OK!");
				  subSeqNameList=null;
				}
				seqObjList=null;
				seqNameList=null;
		    }
		}catch(Exception e){
		  System.out.println(e);
		}
		
		return subSeqObjList;
 } 
 
 public static void extractSubSeq(String seqFile,String subSeqNameFile,
		 int seqNameColIdx,String outSeqFile){
	    
	    List<String> totalSubSeqNames=getRowName(subSeqNameFile,seqNameColIdx);
	    int subSeqNum=totalSubSeqNames.size();

	    int total;
	    if(isFASTASeq(seqFile))
	      total=SeqOperation.getFASTASeqNum(seqFile);
	    else if(isFASTQSeq(seqFile))
		  total=SeqOperation.getFASTQSeqNum(seqFile);
		else
		  return ;
		
		int stopPos=0;
	    int step=10000;
	    int start;
	    int end;
	    List<String> subSeqNameList;
	   	try{
		    String seqIdentifier;
		    String seq;
			String seqQuality="";
			BufferedWriter writer=null;
			writer=new BufferedWriter(new FileWriter(outSeqFile));
			
	   		for(int f=0;f<total;f=f+splitStep){
		  		stopPos=Math.min(total, (f+splitStep));
		   		List<SeqInfo> seqObjList = null;
				if(isFASTASeq(seqFile)){
				   seqObjList=getFASTASeqObj(seqFile,f,stopPos);
				}else if(isFASTQSeq(seqFile)){
				   seqObjList=getFASTQSeqObj(seqFile,f,stopPos);
				}
				List<String> seqNameList=new ArrayList<String>();
				for(SeqInfo seqObj:seqObjList){
				  seqNameList.add(seqObj.seqName);
				}	
				
				for(start=0;start<subSeqNum; start=start+step){				
				  end=Math.min(subSeqNum, (start+step));
				  subSeqNameList=getRowName(totalSubSeqNames,start,end);
			
				  if(isFASTQSeq(outSeqFile)){
						int i=0;
						for(String name: subSeqNameList){
						    i=seqNameList.indexOf(name);
							if(i>=0){
								seqIdentifier=seqObjList.get(i).seqIdentifier;
								seq=seqObjList.get(i).seq;
								seqQuality=seqObjList.get(i).seqQualityEncode;
								
								writer.write("@"+seqIdentifier);
								writer.newLine();
								writer.write(seq);
								writer.newLine();
								writer.write("+");
								writer.newLine();
								writer.write(seqQuality);
								writer.newLine();
								writer.flush();
							}
						}
				  }else{ // fasta format
						int i=0;
						for(String name: subSeqNameList){
						    i=seqNameList.indexOf(name);
							if(i>=0){
								seqIdentifier=seqObjList.get(i).seqIdentifier;
								seq=seqObjList.get(i).seq;
								writer.write(">"+seqIdentifier);
								writer.newLine();
								writer.write(seq);
								writer.newLine();
								writer.flush();
							}
						}
				  } 
				  
				  System.out.println((f+1)+"."+(start+1)+"-"+end+" OK!");
				  subSeqNameList=null;
				}
				seqObjList=null;
				seqNameList=null;
		    }
	  	    writer.close();
		}catch(Exception e){
		  System.out.println(e);
		}

 } 
 

 public static void getSubSeq(String seqFile,String subSeqNameFile,int seqNameColIdx,
		 String outSeqFile){
	    
	    String format=FileOperate.getFileFormat(seqFile);
	    if(outSeqFile==null) 
	       outSeqFile=subSeqNameFile.substring(
	        		       0,subSeqNameFile.lastIndexOf(".")
	        		    )+"."+format;
		
	    try{
			List<SeqInfo> seqObjList = null;
			if(isFASTASeq(seqFile)){
			   seqObjList=getFASTASeqObj(seqFile);
			}else if(isFASTQSeq(seqFile)){
			   seqObjList=getFASTQSeqObj(seqFile);
			}
			List<String> subSeqNameList=getRowName(subSeqNameFile,seqNameColIdx);
			getSubSeq(seqObjList,subSeqNameList,outSeqFile);
			
		}catch(Exception e){
		  System.out.println(e);
		}
	
 } 
 
 public static void getSubSeq(List<SeqInfo> seqObjList,List<String> subSeqNameList,
		 String outSeqFile){	
	 
	 try{		
			List<String> seqNameList=new ArrayList<String>();
			for(SeqInfo seqObj:seqObjList){
			  seqNameList.add(seqObj.seqName);
			}		
			
			String seqIdentifier="";
			String seq="";
			String seqQuality="";
			BufferedWriter writer=null;
			writer=new BufferedWriter(new FileWriter(outSeqFile));
			if(isFASTQSeq(outSeqFile)){
				int i=0;
				for(String name: subSeqNameList){
				    i=seqNameList.indexOf(name);
					if(i>=0){
						seqIdentifier=seqObjList.get(i).seqIdentifier;
						seq=seqObjList.get(i).seq;
						seqQuality=seqObjList.get(i).seqQualityEncode;
						
						writer.write("@"+seqIdentifier);
						writer.newLine();
						writer.write(seq);
						writer.newLine();
						writer.write("+");
						writer.newLine();
						writer.write(seqQuality);
						writer.newLine();
						writer.flush();
					}
				}
			}else{ // fasta format
				int i=0;
				for(String name: subSeqNameList){
				    i=seqNameList.indexOf(name);
					if(i>=0){
						seqIdentifier=seqObjList.get(i).seqIdentifier;
						seq=seqObjList.get(i).seq;
						writer.write(">"+seqIdentifier);
						writer.newLine();
						writer.write(seq);
						writer.newLine();
						writer.flush();
					}
				}
			} 
			seqObjList=null;
			seqNameList=null;
			subSeqNameList=null;
			writer.close();
			writer=null;
			
     }catch(Exception e){
		  System.out.println(e);
	 }
	
 }
 
 public static void saveSeqList(List<SeqInfo> seqObjList, String outSeqFile){
	 
	 if(isFASTAFormat(outSeqFile)){		
		 saveSeqInfoAsFASTA(seqObjList,outSeqFile);		 
	 }else if(isFASTQFormat(outSeqFile)){		 
		 saveSeqInfoAsFASTQ(seqObjList,outSeqFile);		 
	 }		
 }
 
 public static void saveSeqInfoAsFASTA(List<SeqInfo> seqObjList, String outSeqFile){
	 
	 try{
		String seqIdentifier="";
		String seq="";
		BufferedWriter writer=null;
		if(seqObjList!=null && outSeqFile!=null){
		  writer=new BufferedWriter(new FileWriter(outSeqFile));
			
		  for(SeqInfo seqInfo : seqObjList){		
			seqIdentifier=seqInfo.seqIdentifier;
			seq=seqInfo.seq;
			writer.write(">"+seqIdentifier);
			writer.newLine();
			//writer.flush(); 
			writer.write(seq);
			writer.newLine();
			writer.flush();				
		  }

		  writer.close();
		  writer=null;
		}else{
		  System.err.println("null SeqInfo or outSeqFile!"); 
		}	
	}catch(Exception e){
		  System.out.println(e);
	}	 
 }
 
 public static void saveSeqInfoAsFASTA(SeqInfo seqInfo, String outSeqFile){
	 
	 try{
		String seqName="";
		String seq="";
		BufferedWriter writer=null;
		if(seqInfo!=null && outSeqFile!=null){
		    writer=new BufferedWriter(new FileWriter(outSeqFile));
		
		    seqName=seqInfo.seqName;
			seq=seqInfo.seq;
			writer.write(">"+seqName+" "+seq.length());
			writer.newLine();
			//writer.flush(); 
			writer.write(seq);
			writer.newLine();
			writer.flush();				
		
		    writer.close();
		    writer=null;
	   }else{
		  System.err.println("Empty SeqInfo or outSeqFile!"); 
	   }	
	}catch(Exception e){
		  System.out.println(e);
	}	 
 }
 
 public static void saveSeqInfoAsFASTQ(List<SeqInfo> seqObjList, String outSeqFile){
	 
	 try{
		String seqIdentifier="";
		String seq="";
		String seqQuality="";
		BufferedWriter writer=null;
		if(seqObjList!=null && outSeqFile!=null){
		  writer=new BufferedWriter(new FileWriter(outSeqFile));
			
		  for(SeqInfo seqInfo : seqObjList){		
			    seqIdentifier=seqInfo.seqIdentifier;
				seq=seqInfo.seq;
				seqQuality=seqInfo.seqQualityEncode;
				
				writer.write("@"+seqIdentifier);
				writer.newLine();
				writer.write(seq);
				writer.newLine();
				writer.write("+");
				writer.newLine();
				writer.write(seqQuality);
				writer.newLine();
				writer.flush();			
		  }

		  writer.close();
		  writer=null;
		}else{
		  System.err.println("null SeqInfo or outSeqFile!"); 
		}	
	}catch(Exception e){
		  System.out.println(e);
	}	 
 }
 
 public static void saveSeqInfoAsFASTQ(SeqInfo seqInfo, String outSeqFile){
	 
	 try{
		String seqName="";
		String seq="";
		String seqQuality="";
		BufferedWriter writer=null;
		if(seqInfo!=null && outSeqFile!=null){
		    writer=new BufferedWriter(new FileWriter(outSeqFile));
		
		    seqName=seqInfo.seqName;
			seq=seqInfo.seq;
			seqQuality=seqInfo.seqQualityEncode;
			
			writer.write("@"+seqName+" "+seq.length());
			writer.newLine();
			writer.write(seq);
			writer.newLine();
			writer.write("+");
			writer.newLine();
			writer.write(seqQuality);
			writer.newLine();
			writer.flush();	
			
			writer.close();
			writer=null;
	   }else{
		  System.err.println("Empty SeqInfo or outSeqFile!"); 
	   }	
	}catch(Exception e){
		  System.out.println(e);
	}	 
 }
  
  
 public static void saveSeqAsFASTA(String seq,String seqName,String outSeqFile){
	 
	 try{
	
		BufferedWriter writer=null;
		if(seq!=null && outSeqFile!=null){
		    writer=new BufferedWriter(new FileWriter(outSeqFile));		
            if(seqName==null) seqName="NA";
			writer.write(">"+seqName+" "+seq.length());
			writer.newLine();
			//writer.flush(); 
			writer.write(seq);
			writer.newLine();
			writer.flush();				
		
		    writer.close();
		    writer=null;
	   }else{
		  System.err.println("Empty SeqInfo or outSeqFile!"); 
	   }	
	}catch(Exception e){
		  System.out.println(e);
	}	 
 }
 
 public static List<SeqInfo> combineSeqList( List<SeqInfo> seqList1, 
		 List<SeqInfo> seqList2){
	 
	 List<SeqInfo> seqObjList = new ArrayList<SeqInfo>();
	 for(SeqInfo seqInfo:seqList1){
		 seqObjList.add(seqInfo);
	 }
	 
	 for(SeqInfo seqInfo:seqList2){
		 seqObjList.add(seqInfo);
	 }
	 
	 return seqObjList;
 }

 public static List<String> combineStrList( List<String> strList1, 
		 List<String> strList2){
	 
	 List<String> strList = new ArrayList<String>();
	 for(String str:strList1){
		 if(!strList.contains(str))
		    strList.add(str);
	 }
	 
	 for(String str:strList2){
		 if(!strList.contains(str))
			    strList.add(str);
	 }
	 
	 return strList;
 }
    
 public static List<String> loadChrSeqList(String chrSeqFaFile){
		
		List<String> chrSeqList=new ArrayList<String>();

		try{    
	        BufferedReader br;              
	        br = new BufferedReader(new FileReader(chrSeqFaFile));
		    String line;
			line = br.readLine();				
			while(true){           
		       if (line == null) break;
			   if(line.indexOf(">")>=0){			     
				 line=br.readLine();			
				 if (line == null) break;
				 while(line.indexOf(">")<0){
				  chrSeqList.add(line);
				  line=br.readLine();
				  if (line == null) break;
				 }
				
			   }				 
			}           		   
			br.close();
			br=null;
		
		}catch(IOException e){
	        System.out.println(e);
	    }
				
	    return chrSeqList;
		
 }
    
 public static void setChrSiteSeq(List<SeqChrAlignSite> chrAlignSites,
		  int siteDownStreamSize,int siteUpStreamSize,List<String> chrGenomeLineSeqs){
	  
	  for(ChrSite chrSite:chrAlignSites){
		  chrSite.seq=getChrSiteSeq(chrGenomeLineSeqs,chrSite.chrStart,chrSite.chrEnd,
				  chrSite.strand,siteDownStreamSize, siteUpStreamSize);
	  }
 }
  
 public static String getChrSiteSeq(List<String> chrGenomeLineSeqs,int siteStart,int siteEnd,
		  String siteStrand,int siteDownStreamSize,int siteUpStreamSize){

		String siteDownStreamSeq = null;
		String siteUpStreamSeq = null;
		String siteAroundSeq = null;
		
		int coveredDownFaRow;
		int coveredUpFaRow;
		int seqRowStart;
		int seqRowEnd;
		
		int seqSizePerLine=chrGenomeLineSeqs.get(0).length();		
	
		siteDownStreamSize=(siteEnd-siteStart+1)+siteDownStreamSize;
		coveredDownFaRow=siteDownStreamSize/seqSizePerLine+1;
		coveredUpFaRow=siteUpStreamSize/seqSizePerLine+1;	
		
		//System.out.println("GC Content for "+chrName);	  	
		int start=siteStart;
		if(siteStrand.equals("+")) 
		   start=siteStart;
		else if(siteStrand.equals("-")) 
		   start=siteEnd;
				
		double eventFaRowCoord;
		int eventFaRowCoord_iPart;
		double eventFaRowCoord_fPart;
		int eventFaRowCoord_fPart_pos;
		
		eventFaRowCoord=(1.0d*start)/(seqSizePerLine*1.0d);
		eventFaRowCoord_iPart=(int) Math.floor(eventFaRowCoord);
		eventFaRowCoord_fPart=eventFaRowCoord-Math.floor(eventFaRowCoord);
		eventFaRowCoord_fPart_pos=(int) (seqSizePerLine*eventFaRowCoord_fPart);
		
		if(siteStrand.equals("+")){			  
			siteUpStreamSeq="";
			seqRowStart=eventFaRowCoord_iPart-coveredUpFaRow;
			if(seqRowStart<0) seqRowStart=0;
			for(int n=seqRowStart;n<eventFaRowCoord_iPart;n++){
			  siteUpStreamSeq=siteUpStreamSeq+chrGenomeLineSeqs.get(n);
			}
			siteUpStreamSeq=siteUpStreamSeq+chrGenomeLineSeqs.get(
					eventFaRowCoord_iPart).substring(0,eventFaRowCoord_fPart_pos);
			seqRowStart=siteUpStreamSeq.length()-1-siteUpStreamSize;
			if(seqRowStart<0) seqRowStart=0;
			siteUpStreamSeq=siteUpStreamSeq.substring(seqRowStart,siteUpStreamSeq.length());
				
			siteDownStreamSeq="";
			siteDownStreamSeq=chrGenomeLineSeqs.get(eventFaRowCoord_iPart).substring(
					eventFaRowCoord_fPart_pos,seqSizePerLine);
			seqRowEnd=eventFaRowCoord_iPart+1+coveredDownFaRow;
			if(seqRowEnd>chrGenomeLineSeqs.size()) seqRowEnd=chrGenomeLineSeqs.size();
			for(int n=eventFaRowCoord_iPart+1;n<seqRowEnd;n++){
				siteDownStreamSeq=siteDownStreamSeq+chrGenomeLineSeqs.get(n);
			}
			if(siteDownStreamSeq.length()>siteDownStreamSize)
			  siteDownStreamSeq=siteDownStreamSeq.substring(0,siteDownStreamSize);
			else
			  siteDownStreamSeq=siteDownStreamSeq.substring(0,siteDownStreamSeq.length());	
			 
			siteAroundSeq=siteUpStreamSeq+siteDownStreamSeq;
		 
		}else if(siteStrand.equals("-")){				  
			siteDownStreamSeq="";
			seqRowStart=eventFaRowCoord_iPart-coveredDownFaRow;
			if(seqRowStart<0) seqRowStart=0;
			for(int n=seqRowStart;n<eventFaRowCoord_iPart;n++){
			  siteDownStreamSeq=siteDownStreamSeq+chrGenomeLineSeqs.get(n);
			}
			siteDownStreamSeq=siteDownStreamSeq+chrGenomeLineSeqs.get(
					eventFaRowCoord_iPart).substring(0,eventFaRowCoord_fPart_pos);
			seqRowStart=siteDownStreamSeq.length()-1-siteDownStreamSize;
			if(seqRowStart<0) seqRowStart=0;
			siteDownStreamSeq=siteDownStreamSeq.substring(
					seqRowStart,siteDownStreamSeq.length());
					
			siteUpStreamSeq="";
			siteUpStreamSeq=chrGenomeLineSeqs.get(eventFaRowCoord_iPart).substring(
					eventFaRowCoord_fPart_pos,seqSizePerLine);
			seqRowEnd=eventFaRowCoord_iPart+1+coveredUpFaRow;
			if(seqRowEnd>chrGenomeLineSeqs.size()) seqRowEnd=chrGenomeLineSeqs.size();
			for(int n=eventFaRowCoord_iPart+1;n<seqRowEnd;n++){
			  siteUpStreamSeq=siteUpStreamSeq+chrGenomeLineSeqs.get(n);
			}
			if(siteUpStreamSeq.length()>siteUpStreamSize)
			  siteUpStreamSeq=siteUpStreamSeq.substring(0,siteUpStreamSize);
			else
			  siteUpStreamSeq=siteUpStreamSeq.substring(0,siteUpStreamSeq.length());	
								    
			siteAroundSeq=siteDownStreamSeq+siteUpStreamSeq;
			
		}
		
		return siteAroundSeq;

 }
  
 public static List<String> getChrGenomeLineSeqs(String chrSeqFaFile){
		
		List<String> chrGenomeLineSeqs;	
		chrGenomeLineSeqs=loadChrSeqList(chrSeqFaFile);		
		
		return chrGenomeLineSeqs;
		
 }

}