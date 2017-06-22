#! /usr/bin/perl
use Cwd;
use POSIX;
#############################################
$numArgs = $#ARGV +1;
$ARGV[$argnum];

if($ARGV[0] eq "--help" || $ARGV[0] eq ""){
    printf("\nThis code requires one input option. The syntax is:./todo.pl [OPTION]");
    printf("\nPlease choose from the following options:");
    printf("\n\n./todo.pl --help                                   Prints this message");
    printf("\n\n./todo.pl --Submit <python_cfg> <InputPar.txt>     Submit crab jobs to the grid");
    printf("\n                                                   <InputPar.txt> contains input command template.");
    printf("\n                                                   Option: --n number of jobs to submit [njobs=0 by default].");
    printf("\n\nNotes: The <InputPar.dat> format requires 1 input per line where datasetpath = <datasetpath> start a new dataset.");
    printf("\nThe current input options are:");
    printf("\n    DataType = <DataType>");
    printf("\n    datasetpath = <datasetpath>");
    printf("\n    dbs_url = <dbs_url>");
    printf("\n    publish_data_name = <publish_data_name>");
    printf("\n    output_file = <output_file>");
    printf("\n    lumi_mask = <lumi_mask>");
    printf("\n    total_number_of_lumis = <total_number_of_lumis>");
    printf("\n    total_number_of_events = <total_number_of_events>");
    printf("\n    number_of_jobs = <number_of_jobs>");
    printf("\n    se_white_list = <se_white_list>");
    printf("\n    se_black_list = <se_black_list>");
    printf("\nList of DataTypes:\n");
    printf("    data, h_tautau, hpm_taunu, ttbar, w_lnu, w_enu, w_munu, w_taunu, dy_ll, dy_ee, dy_mumu, dy_tautau\n");
    printf("    ZZ, WW, WZ, qcd\n");
    printf("\n\n./todo.pl --SkimSummary <InputPar.txt> <CodeDir> Produces the SkimSummary.log to summarize the skim information."); 
    printf("\n                                                   This is run after retreiving the output from the job. The output");
    printf("\n                                                   goes in Code/InputData/.\n");
    printf("\n\n./todo.pl --CheckOutputLog <InputPar.txt>  Reads the output of jobs received by running getlog."); 
    printf("\n                                                   Compares the EventCounter line between retiries and copies failed jobs to subdir BackupRetries");
    printf("\n                                                   Use this option before running  --SkimSummary.\n");
    printf("\n./todo.pl --CheckandCleanOutput <InputPar.txt>     Checks number fo output files (crab can produce duplicate files when resubmitting)");
    printf("\n                                                   OutputFilesfromDisk.log - List of output Files on Disk.");
    printf("\n                                                   OutputFilesfromlog.log - List of output Files from log files.");
    printf("\n                                                   DiffLogAndDisk.log - List of files which are on Disk but not in log files"); 
    printf("\n\n");
    exit(0); 
} 

######################################
$InputFile=$ARGV[1];
$njobs=0;
for($l=3;$l<$numArgs; $l++){
    if($ARGV[l] eq "--n"){
	$l++;
	$njobs=$ARGV[l];
    }
}

$time= strftime("_%h_%d_%Y_hr%H_min%M",localtime);

if( $ARGV[0] eq "--Submit" ){

    #organize the Lumi_XYZ.root file
    system(sprintf("rm createandsubmittest; touch createandsubmittest"));
    system(sprintf("rm submitall; touch submitall"));
    system(sprintf("rm getoutput; touch getoutput"));    
    system(sprintf("rm forceResubmit; touch forceResubmit"));
    system(sprintf("rm resubmit; touch resubmit"));
    system(sprintf("rm check; touch check"));
    system(sprintf("rm CLEAN; touch CLEAN"));
    $pythonfile=$ARGV[1];
    $TempDataSetFile=$ARGV[2];
    # Open ListofFile.txt
    @datasetpath;
    @dbs_url;
    @publish_data_name;
    @output_file;
    @lumi_mask;
    @total_number_of_lumis;
    @total_number_of_events;
    @number_of_jobs;
    @se_white_list;
    @se_black_list;
    @DataType;
    @globaltag;
    @pileupfile;
    @preselection;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile! [ABORTING]");
    $idx=-1;
    while ($item = <DAT>) {
	chomp($item);
	($a,$b,$c,$d)=split(/ /,$item);
	if($a eq "DataType"){
	    $idx++;
	    push(@datasetpath,"");
	    push(@dbs_url,"");
	    push(@publish_data_name,"");
	    push(@output_file,"");
	    push(@lumi_mask,"");
	    push(@total_number_of_lumis,"");
	    push(@total_number_of_events,"");
	    push(@number_of_jobs,"");
	    push(@se_white_list,"");
	    push(@se_black_list,"");
	    push(@DataType,$c);
	    push(@globaltag,"");
            push(@pileupfile,"");
            push(@preselection,"");
	}
        if($a eq "process.GlobalTag.globaltag"){
            $globaltag[$idx]=$item;
        }
        if($a eq "datasetpath"){
            $datasetpath[$idx]=$item;
        }
	if($a eq "dbs_url"){
	    $dbs_url[$idx]=$item;
	}
        if($a eq "publish_data_name"){
            $publish_data_name[$idx]=$item . $time;
	}
	if($a eq "output_file"){
	    $output_file[$idx]=$item;
	}
	if($a eq "lumi_mask"){
	    $lumi_mask[$idx]=$item;
	}
	if($a eq "total_number_of_lumis"){
	    $total_number_of_lumis[$idx]=$item;
	}
	if($a eq "total_number_of_events"){
	    $total_number_of_events[$idx]=$item;
	}
	if($a eq "number_of_jobs"){
	    $number_of_jobs[$idx]=$item;
	}
        if($a eq "se_white_list"){
            $se_white_list[$idx]=$item;
        }
        if($a eq "se_black_list"){
            $se_black_list[$idx]=$item;
        }
	if($a eq "Pile_Up_File"){
	    $pileupfile[$idx]=$c;
	}
	if($a eq "preselection"){
	    $preselection[$idx]=$c;
	}
    }
    close(DAT);
    ## create crab files and submit 
    $idx=0;
    foreach $data (@DataType){
	$dir=$DataType[$idx];
	$dir=~ s/.root/_CRAB/g;
	$dir=~ s/DataType =/ /g;
	$dir.=sprintf("%d", $idx);
	printf("\ncreating dir: $dir\n");
	system(sprintf("rm $dir -rf; mkdir $dir; cp crab_TEMPLATE.cfg  $dir/crab.cfg;cp $pythonfile $dir/HLT_Tau_Ntuple_cfg.py"));
	system(sprintf("cp ../../data/$pileupfile[$idx] $dir"));
	system(sprintf("./subs \"<DataType>\"                \"$DataType[$idx]\"                     $dir/HLT_Tau_Ntuple_cfg.py"));
	system(sprintf("./subs \"<globaltag>\"               \"$globaltag[$idx]\"                    $dir/HLT_Tau_Ntuple_cfg.py"));
	system(sprintf("./subs \"<Pile_Up_File>\"            \"$pileupfile[$idx]\"                   $dir/HLT_Tau_Ntuple_cfg.py"));
	system(sprintf("./subs \"<PRESELECTION>\"            \"$preselection[$idx]\"                 $dir/HLT_Tau_Ntuple_cfg.py"));
	system(sprintf("./subs \"<datasetpath>\"             \"$datasetpath[$idx]\"                  $dir/HLT_Tau_Ntuple_cfg.py"));
	system(sprintf("./subs \"<Pile_Up_File>\"            \"$pileupfile[$idx]\"                   $dir/crab.cfg"));
	system(sprintf("./subs \"<datasetpath>\"             \"$datasetpath[$idx]\"                  $dir/crab.cfg"));
	system(sprintf("./subs \"<dbs_url>\"                 \"$dbs_url[$idx] \"                     $dir/crab.cfg"));
	system(sprintf("./subs \"<publish_data_name>\"       \"$publish_data_name[$idx] \"           $dir/crab.cfg"));
	system(sprintf("./subs \"<output_file>\"             \"$output_file[$idx] \"                 $dir/crab.cfg"));
	if($lumi_mask[$idx] ne "none"){
	    $lumifile=$lumi_mask[$idx];
	    $lumifile=~ s/lumi_mask =/ /g;
	    system(sprintf("cp $lumifile $dir"));
	    system(sprintf("./subs \"<lumi_mask>\"              \"$lumi_mask[$idx] \"                $dir/crab.cfg"));
	    system(sprintf("./subs \"<total_number_of_lumis>\"  \"$total_number_of_lumis[$idx] \"    $dir/crab.cfg"));
	}
	system(sprintf("./subs \"<total_number_of_events>\" \"$total_number_of_events[$idx] \"       $dir/crab.cfg"));
	system(sprintf("./subs \"<number_of_jobs>\"         \"$number_of_jobs[$idx] \"               $dir/crab.cfg"));
	system(sprintf("./subs \"<number_of_jobs>\"         \"$number_of_jobs[$idx] \"               $dir/crab.cfg"));
	if($se_white_list[$idx] ne "none" || $se_white_list[$idx] ne ""){
	    system(sprintf("./subs \"<se_white_list>\"          \"$se_white_list[$idx] \"            $dir/crab.cfg"));
	}
	if($se_black_list[$idx] ne "none" || $se_black_list[$idx] ne ""){
            system(sprintf("./subs \"<se_black_list>\"          \"$se_black_list[$idx] \"            $dir/crab.cfg"));
        }
	if($njobs !=0){
	    system(sprintf("cd $dir ; crab -create -submit $njobs ; cd .."));
	}
	else{
	    system(sprintf("echo 'cd $dir ; crab -create ; crab -submit 1; cd ..' >>  createandsubmittest \n"));
	    system(sprintf("echo 'cd $dir ; crab -submit ; cd ..' >>  submitall \n"));
	    system(sprintf("echo 'cd $dir;  crab -status; crab -getoutput; cd ..' >> getoutput \n"));
	    system(sprintf("echo \"cd $dir; crab -status | grep -A 2 \\\"resubmit\\\" | grep \\\"jobs:\\\" | awk '{print \\\"crab -resubmit \\\" \\\$4}' | tee junk; source junk; cd ..\" >>  resubmit \n"));
	    system(sprintf("echo \"cd $dir; crab -status | grep '>>>>>>>>'; cd .. \">> check"));
	    system(sprintf("echo \"cd $dir; crab -status | grep -A 2 \\\"Jobs with Wrapper Exit Code\\\" | grep \\\"jobs:\\\" | awk '{print \\\"crab -resubmit \\\" \\\$4}' | tee junk;crab -status | grep -A 2 \\\"Jobs with Wrapper Exit Code : 0\\\" | grep \\\"jobs:\\\" | awk '{print \\\"crab -resubmit \\\" \\\$4}' |  tee junk1;  diff junk junk1 | awk '{print \\\$2 \\\" \\\" \\\$3 \\\" \\\" \\\$4}'| tee junk3; source junk3; cd .. \" >> forceResubmit"));
	    system(sprintf("echo \"rm  $dir -rf; \">> CLEAN"));
	}
    	$idx++;  
    }
    printf("The Submission Directories have been set up....\n");
    printf("To create and submit a test version of your crab jobs: source createandsubmittest \n");
    printf("To create and submit the crab jobs: source submitall \n");
    printf("To get the log files from the crab jobs: source getoutput \n");
    printf("NOTE: PLEASE VALIDATE YOUR CODE RUNS BEFORE SUBMITTING ALL JOBS (ie source createandsubmittest)\n");

}



if( $ARGV[0] eq "--SkimSummary" ){
    $TempDataSetFile=$ARGV[1];
    $CodeDir=$ARGV[2];
    @ID;
    @NEvents;
    @NEventsErr;
    @NEvents_sel;
    @NEventsErr_sel;
    @NEvents_noweight;
    @NEvents_noweight_sel;
    @NFiles;
    
    @DataType;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile! [ABORTING]");
    $idx=-1;
    system(sprintf("touch junk0"));
    system(sprintf("touch junk1"));
    while ($item = <DAT>) {
        chomp($item);
	push(@DataType,$item);
	system(sprintf(" ls -l  $item  >>  junk0 "));
	opendir(DIR,"$item/");
	my @dir = readdir DIR;
	foreach my $subitem (@dir) {
	    @dirs = grep {(/crab_/)} readdir(DIR);
	}
	closedir DIR;
#        ($a,$b,$c,$d)=split(/ /,$item);
#        if($a eq "DataType"){
#	    push(@DataType,$c);
#	}
    }

    system(sprintf("rm junk0; "));
    system(sprintf("rm junk1;"));
    system(sprintf("cat junk0 | awk '{print \$9}' >& junk1"));
    close(DAT);
    $idx=0;
    $diridx=0;
    foreach $data (@DataType){
	printf("Looking for: $data \n");
        $datadir=$data;
        $datadir=~ s/.root/_CRAB/g;
        $datadir=~ s/DataType =/ /g;
        $datadir.=sprintf("%d", $diridx);
#DefaultCrab3Area/crab_DefaultReqName/results/

#	print "------------------data dut $datadir";
	$diridx++;
	$idx++;
	$myDIR=getcwd;
	opendir(DIR,"$myDIR/$data/");
#	printf("Searching: $myDIR/$data/ \n");
	system(sprintf("ls $myDIR/$data/"));
	@dirs = grep {(/crab_/)} readdir(DIR);
	closedir DIR;
	foreach $subdir (@dirs){
	    printf("Opening Dir: $myDIR/$data/$subdir\n");
	    opendir(SUBDIR,"$myDIR/$data/$subdir/results/");
	    @files = grep { /job_out/ } readdir(SUBDIR);

	   # print "===================== files  @files\n";
	    foreach $file (@files){
		open(INPUT,"$myDIR/$data/$subdir/results/$file")  || die "can't open log file $myDIR/$data/$subdir/results/$file" ;
		while (<INPUT>) {
		    ($i1,$i2,$i3,$i4,$i5,$i6,$i7,$i8,$i9,$i10,$i11,$i12,$i13,$i14,$i15)=split(/ /,$_);
		    if($i3 eq "[EventCounter-AllEvents]:" || $i3 eq "[EventCounter-BeforeTauNtuple]:"){
			$flag=0;
		      IDLOOP: {
			  foreach $currentid (@ID){
			      if($currentid eq $i4){
				  last IDLOOP;
			      }
			      $flag++;
			  }
		      }
			$size=@ID;
			if($flag==$size){
			    push(@ID,$i4);
			    push(@NEvents,0);
			    push(@NEventsErr,0);
			    push(@NEvents_sel,0);
			    push(@NEventsErr_sel,0);
			    push(@NEvents_noweight,0);
			    push(@NEvents_noweight_sel,0);
			    push(@NFiles,0);
			}
			if($i3 eq "[EventCounter-AllEvents]:"){
			    $NEvents[$flag]+=$i9;
			 #//   print "===================== $i9  $i9\n";
			    $NEvents_noweight[$flag]+=$i15;
			    $tmp=$NEventsErr[$flag];
			    $NEventsErr[$flag]=sqrt($tmp*$tmp+$i11*$i11);
			    $NFiles[$flag]+=1;
			}#== CMSSW: [EventCounter-AllEvents]: 10330533 no of weighted events: 118 +/- 10.8628 no. of events: 118 +/- 10.8628
			if($i1 eq "[EventCounter-BeforeTauNtuple]:"){
			    $NEvents_sel[$flag]+=$i9;
			    $NEvents_noweight_sel[$flag]+=$i15;
			    $tmp=$NEventsErr_sel[$flag];
                            $NEventsErr_sel[$flag]=sqrt($tmp*$tmp+$i11*$i11);
			}
		    }
		}
		close(OUTPUT) ;
		$idx=0;
		system(sprintf("rm SkimSummary.log"));
		foreach $currentid (@ID){
#		    system(sprintf("echo \"ID= %d AllEvt= %f AllEvtErr= %f SelEvt= %f  SelEvtErr= %f AllEvtnoweight= %f SelEvtnoweight= %f Eff(weight)= %f Eff(noweight)= %f Nfile= %d \" >> SkimSummary.log",$ID[$idx],$NEvents[$idx],$NEventsErr[$idx],$NEvents_sel[$idx],$NEventsErr_sel[$idx],$NEvents_noweight[$idx],$NEvents_noweight_sel[$idx], $NEvents_sel[$idx]/$NEvents[$idx],$NEvents_noweight_sel[$idx]/$NEvents_noweight[$idx],$NFiles[$idx]));

		    system(sprintf("echo \"ID= %d AllEvt= %f AllEvtErr= %f SelEvt= %f  SelEvtErr= %f AllEvtnoweight= %f SelEvtnoweight= %f Nfile= %d \" >> SkimSummary.log",$ID[$idx],$NEvents[$idx],$NEventsErr[$idx],$NEvents_sel[$idx],$NEventsErr_sel[$idx],$NEvents_noweight[$idx],$NEvents_noweight_sel[$idx], $NFiles[$idx]));
		    $idx++;
		}
	    }
	}
    }
#    system(sprintf("cp SkimSummary.log $CodeDir/InputData/SkimSummary.log "));
#    printf("SkimSummary.log has been generated. It has been copied to $CodeDir/InputData/SkimSummary.log\n ");
}

if( $ARGV[0] eq "--CheckOutputLog" ){

    $TempDataSetFile=$ARGV[1];
    $CodeDir=$ARGV[2];
    @ID;
    @NEvents;
    @NEventsErr;
    @NEvents_sel;
    @NEventsErr_sel;
    @NEvents_noweight;
    @NEvents_noweight_sel;
    @NFiles;
    
    @DataType;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile! [ABORTING]");
    $idx=-1;
    system(sprintf("touch junk0"));
    system(sprintf("touch junk1"));
    system(sprintf("rm JobOutputsRetriesSummary.log"));
    system(sprintf("touch JobOutputsRetriesSummary.log"));

  
    while ($item = <DAT>) {
        chomp($item);
	push(@DataType,$item);
	system(sprintf(" ls -l  $item  >>  junk0 "));
	opendir(DIR,"$item/");
	my @dir = readdir DIR;
	foreach my $subitem (@dir) {
	    @dirs = grep {(/crab_/)} readdir(DIR);
	}
	closedir DIR;
#        ($a,$b,$c,$d)=split(/ /,$item);
#        if($a eq "DataType"){
#	    push(@DataType,$c);
#	}
    }

    system(sprintf("rm junk0; "));
    system(sprintf("rm junk1;"));
    system(sprintf("cat junk0 | awk '{print \$9}' >& junk1"));
    close(DAT);
    $idx=0;
    $diridx=0;
    foreach $data (@DataType){
	printf("Looking for: $data \n");
        $datadir=$data;
        $datadir=~ s/.root/_CRAB/g;
        $datadir=~ s/DataType =/ /g;
        $datadir.=sprintf("%d", $diridx);
#DefaultCrab3Area/crab_DefaultReqName/results/

#	print "------------------data dut $datadir";
	$diridx++;
	$idx++;
	$myDIR=getcwd;
	opendir(DIR,"$myDIR/$data/");
#	printf("Searching: $myDIR/$data/ \n");
	system(sprintf("ls $myDIR/$data/"));
	@dirs = grep {(/crab_/)} readdir(DIR);
	closedir DIR;
	my @ToMove;
	foreach $subdir (@dirs){
	    printf("Opening Dir: $myDIR/$data/$subdir\n");
	    opendir(SUBDIR,"$myDIR/$data/$subdir/results/");
	    system(sprintf("mkdir $myDIR/$data/$subdir/results/BackupRetries"));
	    system(sprintf("echo \"             \" >> JobOutputsRetriesSummary.log")); 
	    system(sprintf("echo \">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  \" >> JobOutputsRetriesSummary.log")); 
	    system(sprintf("echo \" $data/$subdir/results/  \" >> JobOutputsRetriesSummary.log")); 
	    
	    @files = grep { /job_out/ } readdir(SUBDIR);
	    $maxretry=0;
	    foreach $file (@files){
		($i1,$i2,$i3,$i4)=split(/\./,$file);
		if($i3 > $maxretry) {
		    $maxretry = $i3;
		}
	    }
	    for($l=$maxretry;$l >=0; $l--){
		$k=$l-1;
		opendir(SUBDIR,"$myDIR/$data/$subdir/results/");
		@files_temp = grep { /.$l.txt/ } readdir(SUBDIR);
		opendir(SUBDIR,"$myDIR/$data/$subdir/results/");
		@files_temp_temp = grep { /.$k.txt/ } readdir(SUBDIR);
  		foreach $file_temp (@files_temp){
		    ($it1,$it2,$it3,$it4)=split(/\./,$file_temp);
		    foreach  $file_temp_temp (@files_temp_temp){
			($itt1,$itt2,$itt3,$itt4)=split(/\./,$file_temp_temp);
			if($it2 eq $itt2) 
			{  
			    push @ToMove, $file_temp_temp;
			    system(sprintf("echo \"-------_- $file_temp $file_temp_temp  \" >> JobOutputsRetriesSummary.log")); 
			    system(sprintf("diff $myDIR/$data/$subdir/results/$file_temp $myDIR/$data/$subdir/results/$file_temp_temp | grep \"EventCounter\" >> JobOutputsRetriesSummary.log"));
			    system(sprintf("echo \"      \" >> JobOutputsRetriesSummary.log")); 
			    system(sprintf("echo \"      \" >> JobOutputsRetriesSummary.log")); 
			 
		#	    system(sprintf("mv $myDIR/$data/$subdir/results/$file_temp_temp  $myDIR/$data/$subdir/results/BackupRetries"));

			}
		    }
		}
	    }
	
	    $a = "JobOutputsRetriesSummary.log";
	    open FILE, ">>$a" or die;  
	    print FILE  "backuped  logs:    @ToMove \n"  ;
	    foreach $filetomove (@ToMove){
		system(sprintf("mv $myDIR/$data/$subdir/results/$filetomove   $myDIR/$data/$subdir/results/BackupRetries"));
	    }
	    
	}
    }
    printf("\n Examine the file  JobOutputsRetriesSummary.log  for retries difference, only the last retry is kept, previous retries before the last are moved to subdir BackupRetries  \n ");
    

}
 





if( $ARGV[0] eq "--CheckandCleanOutput" ){
    $TempDataSetFile=$ARGV[1];
    @DataType;
    open(DAT, $TempDataSetFile) || die("Could not open file $TempDataSetFile! [ABORTING]");
    $idx=-1;
    while ($item = <DAT>) {
        chomp($item);
        ($a,$b,$c,$d)=split(/ /,$item);
        if($a eq "DataType"){
            push(@DataType,$c);
        }
    }
    close(DAT);


    foreach $data (@DataType){
        printf("Looking for: $data \n");
        $datadir=$data;
        $datadir=~ s/.root/_CRAB/g;
        $datadir=~ s/DataType =/ /g;
        $datadir.=sprintf("%d", $diridx);
        $diridx++;
        $idx++;
        $myDIR=getcwd;
        opendir(DIR,"$myDIR/$datadir/");
        printf("Searching: $myDIR/$datadir/ \n");
        system(sprintf("ls $myDIR/$datadir/"));
        @dirs = grep {(/crab_/)} readdir(DIR);
        closedir DIR;
        foreach $subdir (@dirs){
	    system(sprintf("rm  $myDIR/$datadir/OutputFilesfromDisk.log"));
	    system(sprintf("rm $myDIR/$datadir/OutputFilesfromlog.log"));
            system(sprintf("rm $myDIR/$datadir/DiffLogAndDisk.log"));
            system(sprintf("rm FileSummary.log"));
            printf("Opening Dir: $myDIR/$datadir/$subdir\n");
            opendir(SUBDIR,"$myDIR/$datadir/$subdir/res/");
            @files = grep { /stdout/ } readdir(SUBDIR);
	    system(sprintf("rm junk1; touch junk1; touch junk5"));	    
            foreach $file (@files){
		system(sprintf("grep \"LFN:\"  $myDIR/$datadir/$subdir/res/$file  | grep -v \"echo\"  >> junk1"));
	    }
		system(sprintf("cat junk1 | awk '{ split(\$2,a,\"/TauNtuple\"); print \"TauNtuple\"a[2] }' | tee $myDIR/$datadir/OutputFilesfromlog.log"));
		system(sprintf("tail -n 1 junk1 | awk '{ split(\$2,a,\"/TauNtuple\"); print \"uberftp grid-ftp.physik.rwth-aachen.de \\\"cd /pnfs/physik.rwth-aachen.de/cms\" a[1] \" ; ls */ \\\" | tee junk3 \"}' > junk2"));

	    $FilesonDisk="$myDIR/$datadir/OutputFilesfromDisk.log";
	    system(sprintf("echo \"grep root junk3 | awk '{print \\\$8 }'| tee $FilesonDisk   \" >> junk2"));
	    system(sprintf("source junk2;"));

	    open(DAT, $FilesonDisk) || die("Could not open file $FilesonDisk!");
            while ($item = <DAT>) {
                chomp($item);
		system(sprintf("tail -n 1 junk1 | awk '{ split(\$2,a,\"/TauNtuple\"); print \"uberftp grid-ftp.physik.rwth-aachen.de \\\"cd /pnfs/physik.rwth-aachen.de/cms\"  a[1] \" ; rm $item \\\" \"}' >> junk5"));

	    }
	    close(DAT);

	    system(sprintf("echo 'N Files found in in$myDIR/$datadir/OutputFilesfromlog.log: ' >> FileSummary.log "));
	    system(sprintf("cat $myDIR/$datadir/OutputFilesfromlog.log | wc -l >> FileSummary.log"));
	    system(sprintf("echo 'N Files found in in$myDIR/$datadir/OutputFilesfromDisk.log: ' >> FileSummary.log "));
	    system(sprintf("cat $myDIR/$datadir/OutputFilesfromDisk.log | wc -l >> FileSummary.log"));

	    @ListfromDisk=();
	    open(DAT, "$myDIR/$datadir/OutputFilesfromDisk.log") || die("Could not open file $myDIR/$datadir/OutputFilesfromDisk.log! [ABORTING]");
	    $idx=-1;
	    while ($item = <DAT>) {
		chomp($item);
		push(@ListfromDisk,$item);
	    }
	    close(DAT);

	    @ListfromLog=();
	    open(DAT, "$myDIR/$datadir/OutputFilesfromlog.log") || die("Could not open file $myDIR/$datadir/OutputFilesfromlog.log! [ABORTING]");
	    $idx=-1;
	    while ($item = <DAT>) {
		chomp($item);
		push(@ListfromLog,$item);
	    }
	    close(DAT);
	   

	    $DiffFile="$myDIR/$datadir/DiffLogAndDisk.log";
	    system(sprintf("rm $DiffFile; touch  $DiffFile")); 
	    foreach $thedisk (@ListfromDisk){
		$match=0;
		foreach $thelog (@ListfromLog){
		    if($thedisk eq $thelog){
			printf("match");
			$match=1;
		    }
		}
		printf("$thedisk $match\n"); 
		if($match != 1){
		    system(sprintf("echo '$thedisk' >> $DiffFile "));
		}
	    }


	    system(sprintf("rm junk6; touch junk6"));
	    system(sprintf("grep \"LFN:\"  $myDIR/$datadir/$subdir/res/$file  >> junk6"));

	    $CleanFile="$myDIR/$datadir/CleanDisk.dat";
	    $redo="$myDIR/$datadir/Redo";
	    system(sprintf("rm $CleanFile; touch $CleanFile ; rm $redo; touch $redo" ));
	    open(DAT, $DiffFile) || die("Could not open file $DiffFile!");
	    while ($item = <DAT>) {
		chomp($item);
		system(sprintf("grep $item junk5 >> $CleanFile "));
		system(sprintf("grep $item junk6 | awk '{ split(\$1,a,\"/CMSSW_\"); print a[2]}'| awk '{ split(\$1,a,\"/.stdout\"); print \" crab -submit \" a[1]}' | >>$redo")); 
	    }
	    close(DAT);
	    
            system(sprintf("rm junk1; touch junk1"));
            system(sprintf("rm junk2; touch junk2"));
            system(sprintf("rm junk3; touch junk3"));
            system(sprintf("rm junk4; touch junk4"));
	    system(sprintf("rm junk5; touch junk5"));
	}
    }
}
