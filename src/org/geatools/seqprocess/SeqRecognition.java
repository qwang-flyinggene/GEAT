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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.geatools.data.structure.BlastInfo;
import org.geatools.data.structure.SeqCompoAlignInfo;
import org.geatools.data.structure.SeqCompoRecognizer;
import org.geatools.data.structure.SeqRocket;
import org.geatools.data.structure.SeqRocketPair;
import org.geatools.operation.FileOperate;

public class SeqRecognition {
    
	List<SeqRocket> seqRockets;
	List<SeqRocketPair> seqPairRockets;
	boolean isRocketsOK=false;	  

	double extendTerritoryRatio=0.3d;
	  
	static String homeDir=".";
	List<String> tmpFiles;
	String tmpDir;
	String dataDir;
	  
	String regexSeqID = "(\\s+SeqID=)(\\d+)(#\\s*)";
	
	public void setHomeDir(String dir){
	    homeDir=dir;
	}	  
	public void setDataDir(String dir){
	    dataDir=dir;
	}	  
	public void setTmpDir(String dir){
	    tmpDir=dir;
	}	  
	public boolean isSeqRocketsOK(){
	    return isRocketsOK;
	}
	public void setSeqRockets(List<SeqRocket> rockets){
	    seqRockets=rockets;
	    if(rockets!=null && rockets.size()>0) isRocketsOK=true;
		else isRocketsOK=false;;
	}
	public List<SeqRocket>  getSeqRockets(){
	    return seqRockets;
	}  
	public void setSeqPairRockets(List<SeqRocketPair> pairRockets){
	    seqPairRockets=pairRockets;
		if(pairRockets!=null && pairRockets.size()>0) isRocketsOK=true;
		else isRocketsOK=false;;
	}
	public List<SeqRocketPair>  getSeqPairRockets(){
	    return seqPairRockets;
	}
	
	void createFASTASeq(String seq, String seqName, String outSeqFile){
		  
			try{    
		        BufferedWriter writer=null;
		        writer=new BufferedWriter(new FileWriter(outSeqFile));
				String nameLine;
			    //int seqID=0;
			
		        nameLine=">"+seqName+" "+seq.length();
				writer.write(nameLine);
				writer.newLine();
				  //writer.flush(); 
				writer.write(seq);
				writer.newLine();
				writer.flush(); 
				
				writer.close();
				
			}catch(IOException e){
		        System.out.println(e);
		    }  
		 
	}
		   
	boolean checkStruct(){
		  
	    boolean isOK=true;
			
		return isOK;
		   
	}

	List<SeqCompoAlignInfo> getSeqObj(String inSeqFile){
		  
		    List<SeqCompoAlignInfo> seqObjList=null;
			
		    if(SeqOperation.isFASTASeq(inSeqFile)){		
			   seqObjList=getSeqObjFromFASTA(inSeqFile);	   
			}else if(SeqOperation.isFASTQSeq(inSeqFile)){		
			   seqObjList=getSeqObjFromFASTQ(inSeqFile);	   
			}
		    
			return seqObjList;
	}
		  
	List<SeqCompoAlignInfo> getSeqObjFromFASTA(String seqFile){
			  
			    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo>();
			    SeqCompoAlignInfo perSeq=new SeqCompoAlignInfo();
				int seqNum=0;
				try{    
					BufferedReader br;              
			        br = new BufferedReader(new FileReader(seqFile));
				    String line;
					String seqIdentifier;
					String seqLine;
					String seqName;
					String [] itemSplited;
					seqObjList=new ArrayList<SeqCompoAlignInfo>();
					seqNum=0;
					line = br.readLine();
					if (line == null){ 
					   br.close();
					   return seqObjList;
					}
					while(line.length()==0 || line.matches("\\s*")){
					   line = br.readLine();
					   if (line == null){ 
						   br.close();
						   return seqObjList;
					   }
					}
					while(true){           
				       if (line == null) break;
				       line=line.trim();
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
						 	  
						 perSeq=new SeqCompoAlignInfo();
						 perSeq.seqIdentifier=seqIdentifier;
						 perSeq.seqName=seqName;
						 perSeq.seqLength=seqLine.length();
						 perSeq.seq=seqLine;			 
						 seqObjList.add(perSeq);
						 perSeq=null;             
					   }			 
					}//while
			        		  
					br.close();
					
				}catch(IOException e){
			        System.out.println(e);
			    }
			    //System.out.println("invalid reads :"+(seqNum-uniSeqName.size()));
				//uniSeqName=null;
			    return seqObjList; 
	}
			  
	List<SeqCompoAlignInfo> getSeqObjFromFASTQ(String seqFile){
			  
			    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo>();
			    SeqCompoAlignInfo perSeq=new SeqCompoAlignInfo();
				int seqNum=0;
				try{    
					BufferedReader br;              
			        br = new BufferedReader(new FileReader(seqFile));
				    String line;
					String seqIdentifier;	
					String seqLine;
					String seqQualityLine;
					String seqName;
					String [] itemSplited;
					seqObjList=new ArrayList<SeqCompoAlignInfo>();
					seqNum=0;
					outerloop:
					while(true){
			           seqNum=seqNum+1;	
			           
					   //for read header/identifier
			           line = br.readLine();			
				       if (line == null) break;
				       while(line.length()==0 || line.matches("\\s*")){
				    	 line = br.readLine();
				    	 if (line == null) break outerloop;
				       }
					   line=line.trim();
			           if(line.indexOf("@")!=0) {
					      System.out.println("Error in reading fastq line: "+line);
						  break;
					   }			   
			      	   seqIdentifier=line.substring(1,line.length());				 
					   itemSplited=seqIdentifier.split("\\s+");
					   seqName=itemSplited[0].trim();
					   
					   //for read sequence
					   line=br.readLine();			 
					   if (line == null) break;
				       while(line.length()==0 || line.matches("\\s*")){
				    	 line = br.readLine();
				    	 if (line == null) break outerloop;
				       }
					   seqLine=line.trim();			 
					   //seqLine=seqLine.replaceAll("N","n");
					   
					   //for read '+' tag
					   line=br.readLine();			  
					   if (line == null) break;
				       while(line.length()==0 || line.matches("\\s*")){
				    	 line = br.readLine();
				    	 if (line == null) break outerloop;
				       }
					   line=line.trim();
					   if(line.indexOf("+")!=0) {
						  System.out.println("Error in reading fastq line: "+line);
						  break;
					   }		      
			           
					   //for read base quality encode
			           line=br.readLine();			 
					   if (line == null) break;
					   while(line.length()==0 || line.matches("\\s*")){
						  line = br.readLine();
						  if (line == null) break outerloop;
					   }
					   seqQualityLine=line.trim();
			           
			           perSeq=new SeqCompoAlignInfo();
			           perSeq.seqIdentifier=seqIdentifier;
					   perSeq.seqName=seqName;
					   perSeq.seqLength=seqLine.length();
					   perSeq.seq=seqLine;				 
				       perSeq.seqQualityEncode=seqQualityLine;			 
					   seqObjList.add(perSeq);
					   perSeq=null;             		   
			    
					}//while
			        	  
					br.close();
					
					line=null;
					seqIdentifier=null;
					seqLine=null;
					seqName=null;
					itemSplited=null;
					
				}catch(IOException e){
			        System.out.println(e);
			    }
			    //System.out.println("invalid reads :"+(seqNum-uniSeqName.size()));
				//uniSeqName=null;
			    return seqObjList; 
	}
		  
	List<SeqCompoAlignInfo> checkSeq(String inSeqFile){
			    
			    List<SeqCompoAlignInfo> seqObjList = null;

				if(SeqOperation.isFASTASeq(inSeqFile)){			  
					seqObjList=checkFASTASeq(inSeqFile);			
				}else if(SeqOperation.isFASTQSeq(inSeqFile)){		  
					seqObjList=checkFASTQSeq(inSeqFile);			
				}
				
				return seqObjList;
	}
		 
	List<SeqCompoAlignInfo> checkFASTASeq(String seqFile){
		  
		    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo>();
		    SeqCompoAlignInfo perSeq;
			int seqNum=0;

			try{    

				BufferedReader br;              
		        br = new BufferedReader(new FileReader(seqFile));
			    String line;
				String seqIdentifier;
				String seqLine;
				String seqName;
				String [] itemSplited;			
				seqNum=0;
				line = br.readLine();
				if (line == null){ 
				   br.close();
				   return seqObjList;
				}
				while(line.length()==0 || line.matches("\\s*")){
				   line = br.readLine();
				   if (line == null){ 
					   br.close();
					   return seqObjList;
				   }
				}
				
				while(true){           
			       if (line == null) break;
			       line=line.trim();
				   if(line.indexOf(">")==0){
				     seqNum=seqNum+1;
				     seqIdentifier=line.substring(1,line.length());			 
					 itemSplited=seqIdentifier.split("\\s+");
					 seqName=itemSplited[0].trim();
					 seqIdentifier=seqName; ////
					 
					 line=br.readLine();
					 seqLine="";
					 if (line == null) break;
					 while(line.indexOf(">")<0){
				      seqLine = seqLine+line.trim();
					  line=br.readLine();
					  if (line == null) break;
					 }
					 seqLine=seqLine.replaceAll("N","n");
				  
					 perSeq=new SeqCompoAlignInfo();
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
		 
	List<SeqCompoAlignInfo> checkFASTQSeq(String seqFile){
		  
		    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo>();
		    SeqCompoAlignInfo perSeq;
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

				seqNum=0;
				outerloop:	
				while(true){  
		           seqNum=seqNum+1;	
		           
				   //for read header/identifier
		           line = br.readLine();			
			       if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   line=line.trim();
		           if(line.indexOf("@")!=0) {
				      System.out.println("Error in reading fastq line: "+line);
					  break;
				   }			   
		      	   seqIdentifier=line.substring(1,line.length());				 
				   itemSplited=seqIdentifier.split("\\s+");
				   seqName=itemSplited[0].trim();
				   seqIdentifier=seqName; ////
				   
				   //for read sequence
				   line=br.readLine();			 
				   if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   seqLine=line.trim();			 
				   seqLine=seqLine.replaceAll("N","n");
				   
				   //for read '+' tag
				   line=br.readLine();			  
				   if (line == null) break;
			       while(line.length()==0 || line.matches("\\s*")){
			    	 line = br.readLine();
			    	 if (line == null) break outerloop;
			       }
				   line=line.trim();
				   if(line.indexOf("+")!=0) {
					  System.out.println("Error in reading fastq line: "+line);
					  break;
				   }		      
		           
				   //for read base quality encode
		           line=br.readLine();			 
				   if (line == null) break;
				   while(line.length()==0 || line.matches("\\s*")){
					  line = br.readLine();
					  if (line == null) break outerloop;
				   }
				   seqQualityLine=line.trim();

				   perSeq=new SeqCompoAlignInfo();
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
	
	void createRecognizerBLASTTarSeq(List<SeqCompoAlignInfo> seqObjList, SeqRocket seqRocket, 
			  SeqCompoRecognizer recognizer, String outSeqFile){
	  
		try{    
	        int trimSEndIndex=-1;
			int trimLeftShift=0;
			int trimSEnd=0;
			for(int i=0;i<seqRocket.seqRecognizer.size();i++){
			  if(seqRocket.seqRecognizer.get(i).done 
					  && seqRocket.seqRecognizer.get(i).leftShiftForNext ){		 
			   trimSEndIndex=seqRocket.seqRecognizer.get(i).index;
			   trimLeftShift=seqRocket.seqRecognizer.get(i).trimLeftShift;
			  }
			}
			
	        int territoryLen=recognizer.territoryLen;
	        boolean leftSubForBlast=recognizer.leftSubForBlast;    
			
			BufferedWriter writer=null;
	        writer=new BufferedWriter(new FileWriter(outSeqFile));
			String nameLine;
			String seqLine="";
			String trimedSeqLine="";
		    int seqID=0;
			int seqLen=0;
	    	if(trimSEndIndex>=0){
			  if(leftSubForBlast){
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  //nameLine=">"+seqObjList.get(seqID).seqName+" "+seqObjList.get(seqID).seqLength;
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
				  trimSEnd=seqObjList.get(seqID).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;	
				  trimedSeqLine="nnnnnnnnnn";			
				  if(trimSEnd+territoryLen<seqLen)
					trimedSeqLine=seqLine.substring(trimSEnd,trimSEnd+territoryLen);
				  else if(trimSEnd<seqLen)
					trimedSeqLine=seqLine.substring(trimSEnd,seqLen);
					 
				  writer.write(trimedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }
			  }else{
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  //nameLine=">"+seqObjList.get(seqID).seqName+" "+seqObjList.get(seqID).seqLength;
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
				  trimSEnd=seqObjList.get(seqID).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;	
				  trimedSeqLine="nnnnnnnnnn";	
				  if(trimSEnd<seqLen) trimedSeqLine=seqLine.substring(trimSEnd,seqLen);			 
				  
				  writer.write(trimedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }		 
			  }
			}else{
			  if(leftSubForBlast){
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();
					
				  if(territoryLen<seqLen)
					trimedSeqLine=seqLine.substring(0,territoryLen);
				  else 
					trimedSeqLine=seqLine.substring(0,seqLen);
					 
				  writer.write(trimedSeqLine);
				  writer.newLine();
				  writer.flush(); 
			    }
			  }else{
			    for(int i=0;i<seqObjList.size();i++){	
				  seqID=i;	
				  nameLine=">"+seqID;
				  writer.write(nameLine);
				  writer.newLine();
				  //writer.flush(); 
				  seqLine=seqObjList.get(seqID).seq;
				  seqLen=seqLine.length();			 
				  
				  writer.write(seqLine);
				  writer.newLine();
				  writer.flush(); 
			    }		  
			  }
			
			}
			writer.close();
			
		}catch(IOException e){
	        System.out.println(e);
	    }  
	 
	}
	
	void createBLASTTarSeq(List<SeqCompoAlignInfo> seqObj,String outSeqFile){
  		
		try{    

			String seqLine;
			BufferedWriter writer=null;
	        writer=new BufferedWriter(new FileWriter(outSeqFile));	
			for(int i=0;i<seqObj.size();i++){     
		      
				 seqLine=seqObj.get(i).seq;
				 writer.write(">"+i);
				 writer.newLine();
				 writer.write(seqLine);
				 writer.newLine();
				 writer.flush(); 
				 
			}           		   
			
			writer.close();
		}catch(IOException e){
	        System.out.println(e);
	    }

    } 

	  
	void createLeftSubBLASTTarSeq(List<SeqCompoAlignInfo> seqObj,int leftSubSeqLen,String outSeqFile){
		  		
			try{    

				String seqLine;
				//String seqName;
				int subSeqLength=0;
				BufferedWriter writer=null;
		        writer=new BufferedWriter(new FileWriter(outSeqFile));	
				for(int i=0;i<seqObj.size();i++){     
			      
					 seqLine=seqObj.get(i).seq;
					 subSeqLength=leftSubSeqLen;
					 if(leftSubSeqLen>seqLine.length()) subSeqLength=seqLine.length();
		             seqLine=seqLine.substring(0,subSeqLength);
					  
					 writer.write(">"+i);
					 writer.newLine();
					 writer.write(seqLine);
					 writer.newLine();
					 writer.flush(); 
					 
				}           		   
				
				writer.close();
			}catch(IOException e){
		        System.out.println(e);
		    }

	} 

	void initSeqAlignArray(List<SeqCompoAlignInfo> seqObjList, int seqTypeNum){  			

		      ArrayList<Integer> seqAlignSStart;	
		      ArrayList<Integer> seqAlignSEnd;  
		    
		      for(int s=0;s<seqObjList.size();s++){   
			        
					 seqAlignSStart=new ArrayList<Integer>();
					 seqAlignSEnd=new ArrayList<Integer>();
			         for(int i=0;i<seqTypeNum;i++){	           
						seqAlignSStart.add(-1);
						seqAlignSEnd.add(-1);
			         }
			         seqObjList.get(s).seqAlignSStart=seqAlignSStart;
			         seqObjList.get(s).seqAlignSEnd=seqAlignSEnd;
		             seqAlignSStart=null;
					 seqAlignSEnd=null;		
				   		 
			  }	//for	
		     
	}
	
	List<SeqCompoAlignInfo> getRecognizedSeq(List<SeqCompoAlignInfo> seqObjList, 
			int compoIndx){
	 
	    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;		
		int alignSStartIndex=compoIndx;		  
		int alignSEndIndex=alignSStartIndex;	
		
		for(int i=0; i<seqObjList.size(); i++){
		
		   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)>=0 
				   && seqObjList.get(i).seqAlignSEnd.get(alignSEndIndex)>=0){
			   
		     seqObj.add(seqObjList.get(i));
		     
		   }
		}
		
		return seqObj; 
	}
	 
	List<SeqCompoAlignInfo> getNoRecognizedSeq(List<SeqCompoAlignInfo> seqObjList,
			int compoIndx){
	 
	    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;	
		
		int alignSStartIndex=compoIndx;		  
		//int alignSEndIndex=alignSStartIndex;
		
		for(int i=0; i<seqObjList.size(); i++){
		   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)==-1)
		     seqObj.add(seqObjList.get(i));
		
		}
		
		return seqObj;
	 
	} 


	List<SeqCompoAlignInfo> getRecognizedSeq(List<SeqCompoAlignInfo> seqObjList, 
			SeqRocket seqRocket,String seqType){
	 
	    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;		
		int alignSStartIndex=seqRocket.seqTypeInfo.seqTypeName.indexOf(seqType);		  
		int alignSEndIndex=alignSStartIndex;	
		
		for(int i=0; i<seqObjList.size(); i++){
		
		   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)>=0 
				   && seqObjList.get(i).seqAlignSEnd.get(alignSEndIndex)>=0){
			   
		     seqObj.add(seqObjList.get(i));
		     
		   }
		}
		
		return seqObj; 
	}
	 
	List<SeqCompoAlignInfo> getNoRecognizedSeq(List<SeqCompoAlignInfo> seqObjList,
			SeqRocket seqRocket,String seqType){
	 
	    List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> () ;	
		
		int alignSStartIndex=seqRocket.seqTypeInfo.seqTypeName.indexOf(seqType);		  
		//int alignSEndIndex=alignSStartIndex;
		
		for(int i=0; i<seqObjList.size(); i++){
		   if(seqObjList.get(i).seqAlignSStart.get(alignSStartIndex)==-1)
		     seqObj.add(seqObjList.get(i));
		
		}
		
		return seqObj;
	 
	} 
	  
	void setSeqExactAlignInfo(List<SeqCompoAlignInfo> seqObjList,SeqRocket seqRocket,
			 String seqType){
	  
	 	try{
			String seqLine;		
			int leftSStartIndex=0;
	        int leftSStart=0;
	        int leftSEnd=0;		
			int trimSEndIndex=-1;
			int trimLeftShift=0;
			int trimSEnd=0;		
			for(int i=0;i<seqRocket.seqRecognizer.size();i++){
			  if(seqRocket.seqRecognizer.get(i).done 
					  && seqRocket.seqRecognizer.get(i).leftShiftForNext ){		 
			    
				  trimSEndIndex=seqRocket.seqRecognizer.get(i).index;
			      trimLeftShift=seqRocket.seqRecognizer.get(i).trimLeftShift;
			  
			  }
			}
			
			int alignSStartIndex=seqRocket.seqTypeInfo.seqTypeName.indexOf(seqType);		  
			int alignSEndIndex=alignSStartIndex;
			
			SeqCompoRecognizer recognizer=seqRocket.seqRecognizer.get(
					 seqRocket.seqTypeInfo.seqTypeName.indexOf(seqType));	
			String recogniSeq=recognizer.seq;
	        int leftMaxStartIndex=recognizer.maxExactStart-1;
	        if (leftMaxStartIndex<0) leftMaxStartIndex=0;
	        
	        if(trimSEndIndex>=0){
			 for(int i=0;i<seqObjList.size();i++){

				seqLine=seqObjList.get(i).seq;				
				leftSStartIndex=seqLine.indexOf(recogniSeq);
			    trimSEnd=seqObjList.get(i).seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   	
				if(leftSStartIndex>=trimSEnd && leftSStartIndex<=trimSEnd+leftMaxStartIndex){
					leftSStart=leftSStartIndex+1;
					leftSEnd=leftSStartIndex+recogniSeq.length();		
				}else{
					leftSStart=-1;
					leftSEnd=-1;
				}
				
				seqObjList.get(i).seqAlignSStart.set(alignSStartIndex,leftSStart);
				seqObjList.get(i).seqAlignSEnd.set(alignSEndIndex,leftSEnd);
				seqObjList.get(i).seqIndex=i;		
			 }
		   }else{
			 for(int i=0;i<seqObjList.size();i++){

				seqLine=seqObjList.get(i).seq;		
				
				leftSStartIndex=seqLine.indexOf(recogniSeq);		   	
				if(leftSStartIndex>=0 && leftSStartIndex<=leftMaxStartIndex){
				//if(redSStartIndex>=0){	 
					leftSStart=leftSStartIndex+1;
					leftSEnd=leftSStartIndex+recogniSeq.length();		
				}else{
					leftSStart=-1;
					leftSEnd=-1;
				}
				
				seqObjList.get(i).seqAlignSStart.set(alignSStartIndex,leftSStart);
				seqObjList.get(i).seqAlignSEnd.set(alignSEndIndex,leftSEnd);
				seqObjList.get(i).seqIndex=i;		
			 }
		   }		
	    }catch(Exception e){
	        System.out.println(e);
	    }

	}
	 

    List<SeqCompoAlignInfo> getLeftSideBLASTSeq(List<SeqCompoAlignInfo> seqObjList, 
   		 SeqCompoRecognizer recognizer, String seqType, SeqRocket seqRocket, 
   		 String blastOutFile){
   	
   		List<SeqCompoAlignInfo> seqObj=new ArrayList<SeqCompoAlignInfo> ();
   		SeqCompoAlignInfo seq;
   		
   		int alignLen=0;
   		int mismatchNum=0;
   		int gapNum=0;
   		int qStart=0;
   		int sStart=0;
   		int sEnd=0;
   		//int sEnd0=0;
   		
   		int minAlignLen=recognizer.minAlignLen;
   		int maxMismatchNum=recognizer.maxMismatchNum;
   		int maxGapNum=recognizer.maxGapNum;
   		int maxQStart=recognizer.maxQStart;
   		int maxSStart=recognizer.maxSStart;	
   		double maxMismatchRatio=recognizer.maxMismatchRatio;
   		double maxGapRatio=recognizer.maxGapRatio;
   		
   		//ArrayList<ArrayList <String>> blastOutList=new ArrayList<ArrayList <String>> ();
   	    List<ArrayList <String>> blastOut=FileOperate.getMatrixFromFile(blastOutFile);
   		List <Integer> seqID=new ArrayList <Integer>();
   		List <Integer> sStartList=new ArrayList <Integer>();
   		List <Integer> sEndList=new ArrayList <Integer>();
   		String readName="";
   		int readID0=-1;
   		int readID=0;
   		
   		for(int i=0; i<blastOut.size();i++){
   		  readName=blastOut.get(i).get(BlastInfo.colSName).trim();
   		  readID=Integer.parseInt(readName);
   		  if(readID>readID0){
   				
   			 alignLen=Integer.parseInt(blastOut.get(i).get(BlastInfo.colAlignLen).trim());
   			 mismatchNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colMismatchNum).trim());
   			 gapNum=Integer.parseInt(blastOut.get(i).get(BlastInfo.colGapNum).trim());
   			 qStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colQStart).trim());
   			 sStart=Integer.parseInt(blastOut.get(i).get(BlastInfo.colSStart).trim());
   			 sEnd=Integer.parseInt(blastOut.get(i).get(BlastInfo.colSEnd).trim());	
   	         /*
   			 if(sStart>sEnd){
   			   sEnd0=sStart;
   			   sStart=sEnd;
   			   sEnd=sEnd0;
   			 }
   			 */
   	         maxMismatchNum=(int) Math.ceil(alignLen*maxMismatchRatio);
   			 maxGapNum=(int) Math.ceil(alignLen*maxGapRatio);	 			
   				
   			 if(alignLen>=minAlignLen && mismatchNum<=maxMismatchNum && gapNum<=maxGapNum 
   						&& qStart<=maxQStart && sStart<=maxSStart && sEnd>=sStart){
   				seqID.add(readID); 
   				sStartList.add(sStart);
   				sEndList.add(sEnd);
   				readID0=readID;
   			    //blastOutList.add(blastOut.get(i));
   			}
   		  }
   		}
   	
   		blastOut=null;		
   		//blastOutList=null;
   		  
   		int alignSStartIndex=seqRocket.seqTypeInfo.seqTypeName.indexOf(seqType);	
   		int alignSEndIndex=alignSStartIndex;	    
   		
   		int trimSEndIndex=-1;
   		int trimLeftShift=0;
   		int trimSEnd=0;
   		for(int i=0;i<seqRocket.seqRecognizer.size();i++){
   		  if(seqRocket.seqRecognizer.get(i).done 
   				  && seqRocket.seqRecognizer.get(i).leftShiftForNext ){		 
   			
   			  trimSEndIndex=seqRocket.seqRecognizer.get(i).index;
   			  trimLeftShift=seqRocket.seqRecognizer.get(i).trimLeftShift;
   		  
   		  }
   		}
   			
   		if(trimSEndIndex>=0){
   			for(int i=0;i<seqID.size();i++){
   			   seq =new SeqCompoAlignInfo();	
   			   seq=seqObjList.get(seqID.get(i));
   			   trimSEnd=seq.seqAlignSEnd.get(trimSEndIndex)-trimLeftShift;		   		   
   			   seq.seqAlignSStart.set(alignSStartIndex,trimSEnd+sStartList.get(i));
   			   seq.seqAlignSEnd.set(alignSEndIndex,trimSEnd+sEndList.get(i));
   			   seqObj.add(seq);
   			   seq=null;
   			}
   		}else{
   			for(int i=0;i<seqID.size();i++){
   			   seq =new SeqCompoAlignInfo();	
   			   seq=seqObjList.get(seqID.get(i));
   			   seq.seqAlignSStart.set(alignSStartIndex,sStartList.get(i));
   			   seq.seqAlignSEnd.set(alignSEndIndex,sEndList.get(i));
   			   seqObj.add(seq);
   			   seq=null;
   			}	
   		}

   		seqID=null;
   		sStartList=null;
   		sEndList=null;
   		
   		return seqObj;
   	  
    } 

	
	List<SeqCompoAlignInfo> combineSeqObj(List<SeqCompoAlignInfo> seqObjList1, 
			 List<SeqCompoAlignInfo> seqObjList2){
	 
	    List<SeqCompoAlignInfo> seqObjList=new ArrayList<SeqCompoAlignInfo> ();
		
		for(int i=0; i<seqObjList1.size();i++){
		   seqObjList.add(seqObjList1.get(i));	   
		}
	    
	    for(int i=0; i<seqObjList2.size();i++){
		   seqObjList.add(seqObjList2.get(i));	   
		}
	    
		return seqObjList;
	}
	 
	public static void setSeqEncodeVal(SeqCompoAlignInfo seqInfo){
		 Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 int baseLen=Math.min(seqInfo.seq.length(), SeqOperation.maxEncodeBaseLen);
		 for (int b = 0; b < baseLen; b++) {
			 ch = seqInfo.seq.charAt(b);			   
			 val =SeqOperation.getBase2NumMap().get(ch);
			 if(val!=null){
				encodeValStr=encodeValStr+val;
			 }
		 }
		 encodeVal=Long.parseLong(encodeValStr);
		 seqInfo.seqNumEncode=encodeVal;
		
	 }
	 
	 public static void setSeqEncodeVal(List<SeqCompoAlignInfo> seqInfoList){
	     for(SeqCompoAlignInfo seqInfo:seqInfoList){
	    	 setSeqEncodeVal(seqInfo);
	     }	
	 }
	 
	 public static void setSeqEncodeVal(SeqCompoAlignInfo seqInfo,int encodeBaseLen){
			 
		     Integer val;
			 char ch;
			 long encodeVal=0;
			 String encodeValStr="";
			 int baseLen=Math.min(seqInfo.seq.length(), encodeBaseLen);
			 for (int b = 0; b < baseLen; b++) {
				 ch = seqInfo.seq.charAt(b);			   
				 val =SeqOperation.getBase2NumMap().get(ch);
				 if(val!=null){
					encodeValStr=encodeValStr+val;
				 }
			 }
			 encodeVal=Long.parseLong(encodeValStr);
			 seqInfo.seqNumEncode=encodeVal;
			
	 }
		 
	 public static void setSeqEncodeVal(List<SeqCompoAlignInfo> seqInfoList, int encodeBaseLen){
		     for(SeqCompoAlignInfo seqInfo:seqInfoList){
		    	 setSeqEncodeVal(seqInfo,encodeBaseLen);
		     }	
	 }

	 public static long getSeqEncodeVal(String seq,int encodeBaseLen){
		 
	     Integer val;
		 char ch;
		 long encodeVal=0;
		 String encodeValStr="";
		 int baseLen=Math.min(seq.length(), encodeBaseLen);
		 for (int b = 0; b < baseLen; b++) {
			 ch = seq.charAt(b);			   
			 val =SeqOperation.getBase2NumMap().get(ch);
			 if(val!=null){
				encodeValStr=encodeValStr+val;
			 }
		 }
		 encodeVal=Long.parseLong(encodeValStr);
		
		 return encodeVal;
		
    }
 
	 
}
