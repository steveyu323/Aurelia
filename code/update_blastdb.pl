<!DOCTYPE html>
<html lang="en">    
<head>
    <meta charset="UTF-8">
	<title>/c++/src/app/blast/update_blastdb.pl</title>
 	<base href="//www.ncbi.nlm.nih.gov/IEB/ToolBox/CPP_DOC/lxr/">
    <meta name="description" content="NCBI C++ Toolkit Cross Reference">
	<meta http-equiv="cache-control" content="max-age=43200"> 
    <link rel="stylesheet" href="//www.ncbi.nlm.nih.gov/corehtml/ncbi.css">
	<link rel="stylesheet" type="text/css" href="lxr.css">
    <script type="text/javascript" src="lxr.js"></script>
</head>
<body>
<header>
    <div class="nav-menu">
        <a href="/">NCBI Home</a><br>
        <a href="/IEB/">IEB Home</a><br>
        <a href="//ncbi.github.io/cxx-toolkit/">C++ Toolkit docs</a><br>
        <a href="/IEB/ToolBox/C_DOC/lxr/source">C Toolkit source browser</a><br>
        <a href="//www.ncbi.nlm.nih.gov/IEB/ToolBox/SB/hbr.html">C Toolkit source browser (2)</a>
    </div>
    <div class="nav-title">
        <h1 class="main">NCBI C++ Toolkit Cross Reference</h1>
        <h3 class="banner"><script type="text/javascript">SVNLocations();</script>&nbsp; <span class="banner"><a class="banner" href="/IEB/ToolBox/CPP_DOC/lxr/source/">c++</a>/&#x200B;<a class="banner" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/">src</a>/&#x200B;<a class="banner" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/">app</a>/&#x200B;<a class="banner" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/">blast</a>/&#x200B;<a class="banner" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl">update_blastdb.pl</a></span></h3>
    </div>
    <div class="nav-mode">
        <div><span class='modes-sel'>Source navigation</span><br><a class='modes' href="/IEB/ToolBox/CPP_DOC/lxr/diff/src/app/blast/update_blastdb.pl">Diff markup</a><br><a class='modes' href="/IEB/ToolBox/CPP_DOC/lxr/ident">Identifier search</a><br><a class="modes" href="/IEB/ToolBox/CPP_DOC/lxr/search">General search</a><br></div>
    </div>
</header>
<hr>
<main>
<div class="filecontent">
<pre class="filecontent-src">
<span class="comment">#!/usr/bin/env perl</span>
<span class="comment"># $Id: update_blastdb.pl 91494 2020-11-04 21:15:59Z camacho $</span>
<span class="comment"># ===========================================================================</span>
<span class="comment">#</span>
<span class="comment">#                            PUBLIC DOMAIN NOTICE</span>
<span class="comment">#               National Center for Biotechnology Information</span>
<span class="comment">#</span>
<span class="comment">#  This software/database is a "United States Government Work" under the</span>
<span class="comment">#  terms of the United States Copyright Act.  It was written as part of</span>
<span class="comment">#  the author's official duties as a United States Government employee and</span>
<span class="comment">#  thus cannot be copyrighted.  This software/database is freely available</span>
<span class="comment">#  to the public for use. The National Library of Medicine and the U.S.</span>
<span class="comment">#  Government have not placed any restriction on its use or reproduction.</span>
<span class="comment">#</span>
<span class="comment">#  Although all reasonable efforts have been taken to ensure the accuracy</span>
<span class="comment">#  and reliability of the software and data, the NLM and the U.S.</span>
<span class="comment">#  Government do not and cannot warrant the performance or results that</span>
<span class="comment">#  may be obtained by using this software or data. The NLM and the U.S.</span>
<span class="comment">#  Government disclaim all warranties, express or implied, including</span>
<span class="comment">#  warranties of performance, merchantability or fitness for any particular</span>
<span class="comment">#  purpose.</span>
<span class="comment">#</span>
<span class="comment">#  Please cite the author in any work or product based on this material.</span>
<span class="comment">#</span>
<span class="comment"># ===========================================================================</span>
<span class="comment">#</span>
<span class="comment"># Author:  Christiam Camacho</span>
<span class="comment">#</span>
<span class="comment"># File Description:</span>
<span class="comment">#   Script to download the pre-formatted BLAST databases.</span>
<span class="comment">#</span>
<span class="comment"># ===========================================================================</span>

<span class='reserved'>use </span>strict;
<span class='reserved'>use </span>warnings;
<span class='reserved'>use </span>Net::FTP;
<span class='reserved'>use </span>Getopt::Long;
<span class='reserved'>use </span>Pod::Usage;
<span class='reserved'>use </span>File::stat;
<span class='reserved'>use </span>Digest::MD5;
<span class='reserved'>use </span>Archive::Tar;
<span class='reserved'>use </span>File::Temp;
<span class='reserved'>use </span>JSON::PP;

<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=NCBI_FTP">NCBI_FTP</a> = &gt; <span class="string">"ftp.ncbi.nlm.nih.gov"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=BLAST_DB_DIR">BLAST_DB_DIR</a> = &gt; <span class="string">"/blast/db"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=USER">USER</a> = &gt; <span class="string">"anonymous"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=PASSWORD">PASSWORD</a> = &gt; <span class="string">"anonymous"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a> = &gt; 0;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=MAX_DOWNLOAD_ATTEMPTS">MAX_DOWNLOAD_ATTEMPTS</a> = &gt; 3;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=EXIT_FAILURE">EXIT_FAILURE</a> = &gt; 1;

<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AWS_URL">AWS_URL</a> = &gt; <span class="string">"http://s3.amazonaws.com"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AMI_URL">AMI_URL</a> = &gt; <span class="string">"http://169.254.169.254/latest/meta-data/local-hostname"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AWS_BUCKET">AWS_BUCKET</a> = &gt; <span class="string">"ncbi-blast-databases"</span>;

<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCS_URL">GCS_URL</a> = &gt; <span class="string">"https://storage.googleapis.com"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCP_URL">GCP_URL</a> = &gt; <span class="string">"http://metadata.google.internal/computeMetadata/v1/instance/id"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCP_BUCKET">GCP_BUCKET</a> = &gt; <span class="string">"blast-db"</span>;

<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=BLASTDB_MANIFEST">BLASTDB_MANIFEST</a> = &gt; <span class="string">"blastdb-manifest.json"</span>;
<span class='reserved'>use </span>constant <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=BLASTDB_MANIFEST_VERSION">BLASTDB_MANIFEST_VERSION</a> = &gt; <span class="string">"1.0"</span>;

<span class="comment"># Process command line options</span>
<span class='reserved'>my</span> $opt_verbose = 1;
<span class='reserved'>my</span> $opt_quiet = 0;
<span class='reserved'>my</span> $opt_force_download = 0;     
<span class='reserved'>my</span> $opt_help = 0;
<span class='reserved'>my</span> $opt_passive = 1;
<span class='reserved'>my</span> $opt_blastdb_ver = <span class='reserved'>undef</span>;
<span class='reserved'>my</span> $opt_timeout = 120;
<span class='reserved'>my</span> $opt_showall = <span class='reserved'>undef</span>;
<span class='reserved'>my</span> $opt_show_version = 0;
<span class='reserved'>my</span> $opt_decompress = 0;
<span class='reserved'>my</span> $opt_source;
<span class='reserved'>my</span> $opt_legacy_exit_code = 0;
<span class='reserved'>my</span> $opt_nt = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_num_cores">get_num_cores</a>();
<span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=result">result</a> = <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GetOptions">GetOptions</a>(<span class="string">"verbose+"</span>          =&gt;  \$opt_verbose,
                        <span class="string">"quiet"</span>             =&gt;  \$opt_quiet,
                        <span class="string">"force"</span>             =&gt;  \$opt_force_download,
                        <span class="string">"passive:s"</span>         =&gt;  \$opt_passive,
                        <span class="string">"timeout=i"</span>         =&gt;  \$opt_timeout,
                        <span class="string">"showall:s"</span>         =&gt;  \$opt_showall,
                        <span class="string">"version"</span>           =&gt;  \$opt_show_version,
                        <span class="string">"blastdb_version:i"</span> =&gt;  \$opt_blastdb_ver,
                        <span class="string">"decompress"</span>        =&gt;  \$opt_decompress,
                        <span class="string">"source=s"</span>          =&gt;  \$opt_source,
                        <span class="string">"num_threads=i"</span>     =&gt;  \$opt_nt,
                        <span class="string">"legacy_exit_code"</span>  =&gt;  \$opt_legacy_exit_code,
                        <span class="string">"help"</span>              =&gt;  \$opt_help);
$opt_verbose = 0 <span class='reserved'>if</span> $opt_quiet;
<span class='reserved'>die</span> <span class="string">"Failed to parse command line options\n"</span> <span class='reserved'>unless</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=result">result</a>;
pod2usage({-exitval =&gt; 0, -verbose =&gt; 2}) <span class='reserved'>if</span> $opt_help;
<span class='reserved'>if</span> (<span class='reserved'>length</span>($opt_passive) <span class='reserved'>and</span> ($opt_passive !~ /1|<span class='reserved'>no</span>/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>)) {
    pod2usage({-exitval =&gt; 1, -verbose =&gt; 0,
            -msg =&gt; <span class="string">"Invalid value for passive option: '$opt_passive'"</span>});
}
pod2usage({-exitval =&gt; 0, -verbose =&gt; 2}) <span class='reserved'>unless</span> (<span class='reserved'>scalar</span> @ARGV <span class='reserved'>or</span> 
                                                  <span class='reserved'>defined</span>($opt_showall) <span class='reserved'>or</span>
                                                  $opt_show_version);
<span class='reserved'>if</span> (<span class='reserved'>defined</span> $opt_blastdb_ver) {
    pod2usage({-exitval =&gt; 1, -verbose =&gt; 0, 
               -msg =&gt; <span class="string">"Invalid BLAST database version: $opt_blastdb_ver. Supported values: 4, 5"</span>}) 
        <span class='reserved'>unless</span> ($opt_blastdb_ver == 4 <span class='reserved'>or</span> $opt_blastdb_ver == 5);
}
pod2usage({-exitval =&gt; 1, -verbose =&gt; 0, -msg =&gt; <span class="string">"Invalid number of threads"</span>}) 
    <span class='reserved'>if</span> ($opt_nt &lt;= 0);
<span class='reserved'>if</span> (<span class='reserved'>length</span>($opt_passive) <span class='reserved'>and</span> $opt_passive =~ /<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=n">n</a>|<span class='reserved'>no</span>/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
    $opt_passive = 0;
} <span class='reserved'>else</span> {
    $opt_passive = 1;
}
<span class='reserved'>my</span> $exit_code = 0;
$|++;

<span class='reserved'>if</span> ($opt_show_version) {
    <span class='reserved'>my</span> $revision = <span class="string">'$Revision: 91494 $'</span>;
    $revision =~ <span class='reserved'>s</span>/\$Revision: | \$//<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=g">g</a>;
    <span class='reserved'>print</span> <span class="string">"$0 version $revision\n"</span>;
    <span class='reserved'>exit</span>($exit_code);
}
<span class='reserved'>my</span> $curl = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_curl_path">get_curl_path</a>();

<span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> = <span class="string">"NCBI"</span>;
<span class="comment"># If provided, the source takes precedence over any attempts to determine the closest location</span>
<span class='reserved'>if</span> (<span class='reserved'>defined</span>($opt_source)) {
    <span class='reserved'>if</span> ($opt_source =~ /^<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ncbi">ncbi</a>/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
        $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> = <span class="string">"NCBI"</span>;
    } <span class='reserved'>elsif</span> ($opt_source =~ /^<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=gc">gc</a>/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
        $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> = <span class="string">"GCP"</span>;
    } <span class='reserved'>elsif</span> ($opt_source =~ /^aws/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
        $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> = <span class="string">"AWS"</span>;
    }
} <span class='reserved'>else</span> {
    <span class="comment"># Try to auto-detect whether we're on the cloud</span>
    <span class='reserved'>if</span> (<span class='reserved'>defined</span>($curl)) {
        <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=tmpfile">tmpfile</a> = <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=File">File</a>::<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=Temp">Temp</a>-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=new">new</a>();
        <span class='reserved'>my</span> $gcp_cmd = <span class="string">"$curl --connect-timeout 3 --retry 3 --retry-max-time 30 -sfo $tmpfile -H 'Metadata-Flavor: Google' "</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCP_URL">GCP_URL</a>;
        <span class='reserved'>my</span> $aws_cmd = <span class="string">"$curl --connect-timeout 3 --retry 3 --retry-max-time 30 -sfo /dev/null "</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AMI_URL">AMI_URL</a>;
        <span class='reserved'>print</span> <span class="string">"$gcp_cmd\n"</span> <span class='reserved'>if</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>;
        <span class='reserved'>if</span> (<span class='reserved'>system</span>($gcp_cmd) == 0) { 
            <span class="comment"># status not always reliable.  Check that curl output is all digits.</span>
            <span class='reserved'>my</span> $tmpfile_content = <span class='reserved'>do</span> { <span class='reserved'>local</span> $/; &lt;$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=tmpfile">tmpfile</a>&gt;};
            <span class='reserved'>print</span> <span class="string">"curl output $tmpfile_content\n"</span> <span class='reserved'>if</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>;
            $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> = <span class="string">"GCP"</span> <span class='reserved'>if</span> ($tmpfile_content =~ <span class="extra">m/^(\d+)$/</span>);
        } <span class='reserved'>elsif</span> (<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>) {
            <span class="comment"># Consult https://ec.haxx.se/usingcurl/usingcurl-returns</span>
            <span class='reserved'>print</span> <span class="string">"curl to GCP metadata server returned "</span>, $?&gt;&gt;8, <span class="string">"\n"</span>;
        }

        <span class='reserved'>print</span> <span class="string">"$aws_cmd\n"</span> <span class='reserved'>if</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>;
        <span class='reserved'>if</span> (<span class='reserved'>system</span>($aws_cmd) == 0) {
            $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> = <span class="string">"AWS"</span>;
        } <span class='reserved'>elsif</span> (<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>) {
            <span class="comment"># Consult https://ec.haxx.se/usingcurl/usingcurl-returns</span>
            <span class='reserved'>print</span> <span class="string">"curl to AWS metadata server returned "</span>, $?&gt;&gt;8, <span class="string">"\n"</span>;
        }
        <span class='reserved'>print</span> <span class="string">"Location is $location\n"</span> <span class='reserved'>if</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>;
    }
}
<span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> =~ /aws|gcp/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a> <span class='reserved'>and</span> <span class='reserved'>not</span> <span class='reserved'>defined</span> $curl) {
    <span class='reserved'>print</span> <span class="string">"Error: $0 depends on curl to fetch data from cloud storage, please install this utility to access these data sources.\n"</span>;
    <span class='reserved'>exit</span>(<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=EXIT_FAILURE">EXIT_FAILURE</a>);
}

<span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>;

<span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> <span class='reserved'>ne</span> <span class="string">"NCBI"</span>) {
    <span class='reserved'>die</span> <span class="string">"Only BLASTDB version 5 is supported at GCP and AWS\n"</span> <span class='reserved'>if</span> (<span class='reserved'>defined</span> $opt_blastdb_ver <span class='reserved'>and</span> $opt_blastdb_ver != 5);
    <span class='reserved'>my</span> $latest_dir = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_latest_dir">get_latest_dir</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a>);
    <span class='reserved'>my</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=json">json</a>, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a>) = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_blastdb_metadata">get_blastdb_metadata</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a>, $latest_dir);
    <span class='reserved'>unless</span> (<span class='reserved'>length</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=json">json</a>)) {
        <span class='reserved'>print</span> STDERR <span class="string">"ERROR: Missing manifest file $url, please report to blast-help\@ncbi.nlm.nih.gov\n"</span>;
        <span class='reserved'>exit</span>(<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=EXIT_FAILURE">EXIT_FAILURE</a>);
    }
    <span class='reserved'>print</span> <span class="string">"Connected to $location\n"</span> <span class='reserved'>if</span> $opt_verbose;
    <span class='reserved'>my</span> $metadata = decode_json($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=json">json</a>);
    <span class='reserved'>unless</span> (<span class='reserved'>exists</span>($$metadata{<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=version">version</a>}) <span class='reserved'>and</span> ($$metadata{<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=version">version</a>} <span class='reserved'>eq</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=BLASTDB_MANIFEST_VERSION">BLASTDB_MANIFEST_VERSION</a>)) {
        <span class='reserved'>print</span> STDERR <span class="string">"ERROR: Invalid version in manifest file $url, please report to blast-help\@ncbi.nlm.nih.gov\n"</span>;
        <span class='reserved'>exit</span>(<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=EXIT_FAILURE">EXIT_FAILURE</a>);
    }
    <span class='reserved'>if</span> (<span class='reserved'>defined</span>($opt_showall)) {
        <span class='reserved'>my</span> $print_header = 1;
        <span class='reserved'>foreach</span> <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a> (<span class='reserved'>sort</span> <span class='reserved'>keys</span> %$metadata) {
            <span class='reserved'>next</span> <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a> =~ /^<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=version">version</a>$/);
            <span class='reserved'>if</span> ($opt_showall =~ /tsv/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
                <span class='reserved'>printf</span>(<span class="string">"%s\t%s\t%9.4f\t%s\n"</span>, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>, $$metadata{$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>}{description}, 
                    $$metadata{$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>}{<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=size">size</a>}, $$metadata{$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>}{last_updated});
            } <span class='reserved'>elsif</span> ($opt_showall =~ /pretty/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
                <span class='reserved'>if</span> ($print_header) {
                    <span class='reserved'>printf</span>(<span class="string">"%-60s %-120s %-11s %15s\n"</span>, <span class="string">"BLASTDB"</span>, 
                        <span class="string">"DESCRIPTION"</span>, <span class="string">"SIZE (GB)"</span>, <span class="string">"LAST_UPDATED"</span>);
                    $print_header = 0;
                }
                <span class='reserved'>printf</span>(<span class="string">"%-60s %-120s %9.4f %15s\n"</span>, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>, $$metadata{$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>}{description}, 
                    $$metadata{$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>}{<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=size">size</a>}, $$metadata{$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>}{last_updated});
            } <span class='reserved'>else</span> {
                <span class='reserved'>print</span> <span class="string">"$db\n"</span>;
            }
        }
    } <span class='reserved'>else</span> {
        <span class='reserved'>my</span> @files2download;
        <span class='reserved'>for</span> <span class='reserved'>my</span> $requested_db (@ARGV) {
            <span class='reserved'>if</span> (<span class='reserved'>exists</span> $$metadata{$requested_db}) {
                <span class='reserved'>push</span> @files2download, @{$$metadata{$requested_db}{<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=files">files</a>}};
            } <span class='reserved'>else</span> {
                <span class='reserved'>print</span> STDERR <span class="string">"Warning: $requested_db does not exist in $location ($latest_dir)\n"</span>;
            }
        }
        <span class='reserved'>if</span> (@files2download) {
            <span class='reserved'>my</span> $gsutil = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_gsutil_path">get_gsutil_path</a>();
            <span class='reserved'>my</span> $awscli = <span class='reserved'>undef</span>; <span class="comment"># &amp;get_awscli_path(); # aws s3 required credentials, fall back to curl</span>
            <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a>;
            <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=fh">fh</a> = <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=File">File</a>::<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=Temp">Temp</a>-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=new">new</a>();
            <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> <span class='reserved'>eq</span> <span class="string">"GCP"</span> <span class='reserved'>and</span> <span class='reserved'>defined</span>($gsutil)) {
                $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> = <span class="string">"$gsutil "</span>;
                <span class='reserved'>if</span> ($opt_nt &gt; 1) {
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">"-m -q cp "</span>;
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">"-o 'GSUtil:parallel_thread_count=1' -o 'GSUtil:parallel_process_count=$opt_nt' "</span>;
                } <span class='reserved'>else</span> {
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">"-q cp "</span>;
                }
                $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class='reserved'>join</span>(<span class="string">" "</span>, @files2download) . <span class="string">" ."</span>;
            } <span class='reserved'>elsif</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> <span class='reserved'>eq</span> <span class="string">"AWS"</span> <span class='reserved'>and</span> <span class='reserved'>defined</span> ($awscli)) {
                <span class='reserved'>my</span> $aws_cmd = <span class="string">"$awscli s3 cp "</span>;
                $aws_cmd .= <span class="string">"--only-show-errors "</span> <span class='reserved'>unless</span> $opt_verbose &gt;= 3;
                <span class='reserved'>print</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=fh">fh</a> <span class='reserved'>join</span>(<span class="string">"\n"</span>, @files2download);
                $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> = <span class="string">"/usr/bin/xargs -P $opt_nt -n 1 -I{}"</span>;
                $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">" -t"</span> <span class='reserved'>if</span> $opt_verbose &gt; 3;
                $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">" $aws_cmd {} ."</span>;
                $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">" &lt;$fh "</span> ;
            } <span class='reserved'>else</span> { <span class="comment"># fall back to  curl</span>
                <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a> = $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=location">location</a> <span class='reserved'>eq</span> <span class="string">"AWS"</span> ? <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AWS_URL">AWS_URL</a> : <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCS_URL">GCS_URL</a>;
                <span class='reserved'>s</span>,gs://,$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a>/, <span class='reserved'>foreach</span> (@files2download);
                <span class='reserved'>s</span>,<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=s3">s3</a>://,$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a>/, <span class='reserved'>foreach</span> (@files2download);
                <span class='reserved'>if</span> ($opt_nt &gt; 1 <span class='reserved'>and</span> <span class='reserved'>-f</span> <span class="string">"/usr/bin/xargs"</span>) {
                    <span class='reserved'>print</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=fh">fh</a> <span class='reserved'>join</span>(<span class="string">"\n"</span>, @files2download);
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> = <span class="string">"/usr/bin/xargs -P $opt_nt -n 1"</span>;
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">" -t"</span> <span class='reserved'>if</span> $opt_verbose &gt; 3;
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">" $curl -sSOR"</span>;
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">" &lt;$fh "</span> ;
                } <span class='reserved'>else</span> {
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> = <span class="string">"$curl -sSR"</span>;
                    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> .= <span class="string">" -O $_"</span> <span class='reserved'>foreach</span> (@files2download);
                }
            }
            <span class='reserved'>print</span> <span class="string">"$cmd\n"</span> <span class='reserved'>if</span> $opt_verbose &gt; 3;
            <span class='reserved'>system</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a>);
        }
    }

} <span class='reserved'>else</span> {
    <span class="comment"># Connect and download files</span>
    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a> = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=connect_to_ftp">connect_to_ftp</a>();
    <span class='reserved'>if</span> (<span class='reserved'>defined</span> $opt_showall) {
        <span class='reserved'>print</span> <span class="string">"$_\n"</span> <span class='reserved'>foreach</span> (<span class='reserved'>sort</span>(&amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_available_databases">get_available_databases</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;ls())));
    } <span class='reserved'>else</span> {
        <span class='reserved'>my</span> @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=files">files</a> = <span class='reserved'>sort</span>(&amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_files_to_download">get_files_to_download</a>());
        <span class='reserved'>my</span> @files2decompress;
        $exit_code = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=download">download</a>(\@<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=files">files</a>, \@files2decompress);
        <span class='reserved'>if</span> ($exit_code == 1) {
            <span class='reserved'>foreach</span> (@files2decompress) {
                $exit_code = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=decompress">decompress</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_">_</a>);
                <span class='reserved'>last</span> <span class='reserved'>if</span> ($exit_code != 1);
            }
        }
        <span class='reserved'>unless</span> ($opt_legacy_exit_code) {
            $exit_code = ($exit_code == 1 ? 0 : $exit_code);
        }
    }
    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;quit();
}

<span class='reserved'>exit</span>($exit_code);

<span class="comment"># Connects to NCBI ftp server</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=connect_to_ftp">connect_to_ftp</a>
{
    <span class='reserved'>my</span> %ftp_opts;
    $ftp_opts{<span class="string">'Passive'</span>} = 1 <span class='reserved'>if</span> $opt_passive;
    $ftp_opts{<span class="string">'Timeout'</span>} = $opt_timeout <span class='reserved'>if</span> ($opt_timeout &gt;= 0);
    $ftp_opts{<span class="string">'Debug'</span>}   = 1 <span class='reserved'>if</span> ($opt_verbose &gt; 1);
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a> = Net::FTP-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=new">new</a>(<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=NCBI_FTP">NCBI_FTP</a>, %ftp_opts)
        <span class='reserved'>or</span> <span class='reserved'>die</span> <span class="string">"Failed to connect to "</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=NCBI_FTP">NCBI_FTP</a> . <span class="string">": $!\n"</span>;
    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=login">login</a>(<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=USER">USER</a>, <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=PASSWORD">PASSWORD</a>) 
        <span class='reserved'>or</span> <span class='reserved'>die</span> <span class="string">"Failed to login to "</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=NCBI_FTP">NCBI_FTP</a> . <span class="string">": $!\n"</span>;
    <span class='reserved'>my</span> $ftp_path = <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=BLAST_DB_DIR">BLAST_DB_DIR</a>;
    $ftp_path .= <span class="string">"/v$opt_blastdb_ver"</span> <span class='reserved'>if</span> (<span class='reserved'>defined</span> $opt_blastdb_ver);
    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;cwd($ftp_path);
    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;binary();
    <span class='reserved'>if</span> ($opt_verbose) {
        <span class='reserved'>if</span> (<span class='reserved'>defined</span> $opt_blastdb_ver) {
            <span class='reserved'>print</span> <span class="string">"Connected to $location; downloading BLASTDBv$opt_blastdb_ver\n"</span>;
        } <span class='reserved'>else</span> {
            <span class='reserved'>print</span> <span class="string">"Connected to $location\n"</span>;
        }
    }
    <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>;
}

<span class="comment"># Gets the list of available databases on NCBI FTP site</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_available_databases">get_available_databases</a>
{
    <span class='reserved'>my</span> @blast_db_files = $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;ls();
    <span class='reserved'>my</span> @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = ();

    <span class='reserved'>foreach</span> (@blast_db_files) {
        <span class='reserved'>next</span> <span class='reserved'>unless</span> (/\.<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=tar">tar</a>\.gz$/);
        <span class='reserved'>push</span> @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>, &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=extract_db_name">extract_db_name</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_">_</a>);
    }
    <span class='reserved'>my</span> %seen = ();
    <span class='reserved'>return</span> <span class='reserved'>grep</span> { ! $seen{$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_">_</a>} ++ } @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>;
}

<span class="comment"># Obtains the list of files to download</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_files_to_download">get_files_to_download</a>
{
    <span class='reserved'>my</span> @blast_db_files = $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;ls();
    <span class='reserved'>my</span> @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = ();

    <span class='reserved'>if</span> ($opt_verbose &gt; 2) {
        <span class='reserved'>print</span> <span class="string">"Found the following files on ftp site:\n"</span>;
        <span class='reserved'>print</span> <span class="string">"$_\n"</span> <span class='reserved'>for</span> (@blast_db_files);
    }

    <span class='reserved'>for</span> <span class='reserved'>my</span> $requested_db (@ARGV) {
        <span class='reserved'>for</span> <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> (@blast_db_files) {
            <span class='reserved'>next</span> <span class='reserved'>unless</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> =~ /\.<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=tar">tar</a>\.gz$/);    
            <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> =~ /^$requested_db\..*/) {
                <span class='reserved'>push</span> @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>;
            }
        }
    }

    <span class='reserved'>if</span> ($opt_verbose) {
        <span class='reserved'>for</span> <span class='reserved'>my</span> $requested_db (@ARGV) {
            <span class='reserved'>unless</span> (<span class='reserved'>grep</span>(/$requested_db/, @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>)) {
                <span class='reserved'>print</span> STDERR <span class="string">"$requested_db not found, skipping.\n"</span> 
            }
        }
    }

    <span class='reserved'>return</span> @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>;
}

<span class="comment"># Download the requested files only if their checksum files are missing or if</span>
<span class="comment"># these (or the archives) are newer in the FTP site. Returns 0 if no files were</span>
<span class="comment"># downloaded, 1 if at least one file was downloaded (so that this can be the</span>
<span class="comment"># application's exit code)</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=download">download</a>($$)
{
    <span class='reserved'>my</span> @requested_dbs = @ARGV;
    <span class='reserved'>my</span> @files2download = @{$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_">_</a>[0]};
    <span class='reserved'>my</span> $files2decompress = $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_">_</a>[1];
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = 0;

    <span class='reserved'>for</span> <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> (@files2download) {

        <span class='reserved'>my</span> $attempts = 0;   <span class="comment"># Download attempts for this file</span>
        <span class='reserved'>if</span> ($opt_verbose <span class='reserved'>and</span> &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=is_multivolume_db">is_multivolume_db</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>) <span class='reserved'>and</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> =~ /\.00\./) {
            <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db_name">db_name</a> = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=extract_db_name">extract_db_name</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>);
            <span class='reserved'>my</span> $nvol = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_num_volumes">get_num_volumes</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db_name">db_name</a>, @files2download);
            <span class='reserved'>print</span> <span class="string">"Downloading $db_name ("</span> . $nvol . <span class="string">" volumes) ...\n"</span> <span class='reserved'>unless</span> ($opt_quiet);
        }

        <span class="comment"># We preserve the checksum files as evidence of the downloaded archive</span>
        <span class='reserved'>my</span> $checksum_file = <span class="string">"$file.md5"</span>;
        <span class='reserved'>my</span> $new_download = (<span class='reserved'>-e</span> $checksum_file ? 0 : 1);
        <span class='reserved'>my</span> $update_available = ($new_download <span class='reserved'>or</span> 
                    ((<span class='reserved'>stat</span>($checksum_file))-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=mtime">mtime</a> &lt; $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;mdtm($checksum_file)));
        <span class='reserved'>if</span> (<span class='reserved'>-e</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> <span class='reserved'>and</span> (<span class='reserved'>stat</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>)-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=mtime">mtime</a> &lt; $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;mdtm($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>))) {
            $update_available = 1;
        }

<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=download_file">download_file</a>:
        <span class='reserved'>if</span> ($opt_force_download <span class='reserved'>or</span> $new_download <span class='reserved'>or</span> $update_available) {
            <span class='reserved'>print</span> <span class="string">"Downloading $file..."</span> <span class='reserved'>if</span> $opt_verbose;
            $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get">get</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>);
            <span class='reserved'>unless</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=ftp">ftp</a>-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get">get</a>($checksum_file)) {
                <span class='reserved'>print</span> STDERR <span class="string">"Failed to download $checksum_file!\n"</span>;
                <span class='reserved'>return</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=EXIT_FAILURE">EXIT_FAILURE</a>;
            }
            <span class='reserved'>my</span> $rmt_digest = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=read_md5_file">read_md5_file</a>($checksum_file);
            <span class='reserved'>my</span> $lcl_digest = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=compute_md5_checksum">compute_md5_checksum</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>);
            <span class='reserved'>print</span> <span class="string">"\nRMT $file Digest $rmt_digest"</span> <span class='reserved'>if</span> (<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>);
            <span class='reserved'>print</span> <span class="string">"\nLCL $file Digest $lcl_digest\n"</span> <span class='reserved'>if</span> (<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>);
            <span class='reserved'>if</span> ($lcl_digest <span class='reserved'>ne</span> $rmt_digest) {
                <span class='reserved'>unlink</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>, $checksum_file;
                <span class='reserved'>if</span> (++$attempts &gt;= <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=MAX_DOWNLOAD_ATTEMPTS">MAX_DOWNLOAD_ATTEMPTS</a>) {
                    <span class='reserved'>print</span> STDERR <span class="string">"too many failures, aborting download!\n"</span>;
                    <span class='reserved'>return</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=EXIT_FAILURE">EXIT_FAILURE</a>;
                } <span class='reserved'>else</span> {
                    <span class='reserved'>print</span> <span class="string">"corrupt download, trying again.\n"</span>;
                    <span class='reserved'>goto</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=download_file">download_file</a>;
                }
            }
            <span class='reserved'>push</span> @$files2decompress, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> <span class='reserved'>if</span> ($opt_decompress);
            <span class='reserved'>print</span> <span class="string">" [OK]\n"</span> <span class='reserved'>if</span> $opt_verbose;
            $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = 1 <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> == 0);
        } <span class='reserved'>else</span> {
            <span class='reserved'>if</span> ($opt_decompress <span class='reserved'>and</span> <span class='reserved'>-f</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>) {
                <span class='reserved'>push</span> @$files2decompress, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>;
                $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = 1;
            } <span class='reserved'>else</span> {
                <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=msg">msg</a> = ($opt_decompress 
                           ? <span class="string">"The contents of $file are up to date in your system."</span> 
                           : <span class="string">"$file is up to date."</span>);
                <span class='reserved'>print</span> <span class="string">"$msg\n"</span> <span class='reserved'>if</span> $opt_verbose;
            }
        }
    }
    <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>;
}

<span class="comment"># Try to decompress using /bin/tar as Archive::Tar is known to be slower (as</span>
<span class="comment"># it's pure perl)</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_decompress_impl">_decompress_impl</a>($)
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>if</span> ($^O <span class='reserved'>eq</span> <span class="string">"cygwin"</span>) {
        <span class='reserved'>local</span> $ENV{PATH} = <span class="string">"/bin:/usr/bin"</span>;
        <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> = <span class="string">"tar -zxf $file 2&gt;/dev/null"</span>;
        <span class='reserved'>return</span> 1 <span class='reserved'>unless</span> (<span class='reserved'>system</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a>));
    }
    <span class='reserved'>unless</span> ($^O =~ /win/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
        <span class='reserved'>local</span> $ENV{PATH} = <span class="string">"/bin:/usr/bin"</span>;
        <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> = <span class="string">"gzip -cd $file 2&gt;/dev/null | tar xf - 2&gt;/dev/null"</span>;
        <span class='reserved'>return</span> 1 <span class='reserved'>unless</span> (<span class='reserved'>system</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a>));
    }
    <span class='reserved'>return</span> Archive::Tar-&gt;extract_archive($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>, 1);
}

<span class="comment"># Decompresses the file passed as its argument</span>
<span class="comment"># Returns 1 on success, and 2 on failure, printing an error to STDERR</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=decompress">decompress</a>($)
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>print</span> <span class="string">"Decompressing $file ..."</span> <span class='reserved'>unless</span> ($opt_quiet);
    <span class='reserved'>my</span> $succeeded = &amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_decompress_impl">_decompress_impl</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>);
    <span class='reserved'>unless</span> ($succeeded) {
        <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=msg">msg</a> = <span class="string">"Failed to decompress $file ($Archive::Tar::error), "</span>;
        $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=msg">msg</a> .= <span class="string">"please do so manually."</span>;
        <span class='reserved'>print</span> STDERR <span class="string">"$msg\n"</span>;
        <span class='reserved'>return</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=EXIT_FAILURE">EXIT_FAILURE</a>;
    }
    <span class='reserved'>unlink</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>;   <span class="comment"># Clean up archive, but preserve the checksum file</span>
    <span class='reserved'>print</span> <span class="string">" [OK]\n"</span> <span class='reserved'>unless</span> ($opt_quiet);
    <span class='reserved'>return</span> 1;
}

<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=compute_md5_checksum">compute_md5_checksum</a>($)
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>my</span> $digest = <span class="string">"N/A"</span>;
    <span class='reserved'>if</span> (<span class='reserved'>open</span>(DOWNLOADED_FILE, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>)) {
        <span class='reserved'>binmode</span>(DOWNLOADED_FILE);
        $digest = Digest::MD5-&gt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=new">new</a>-&gt;addfile(*DOWNLOADED_FILE)-&gt;hexdigest;
        <span class='reserved'>close</span>(DOWNLOADED_FILE);
    }
    <span class='reserved'>return</span> $digest;
}

<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=read_md5_file">read_md5_file</a>($)
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=md5file">md5file</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>open</span>(<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=IN">IN</a>, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=md5file">md5file</a>);
    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_">_</a> = &lt;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=IN">IN</a>&gt;;
    <span class='reserved'>close</span>(<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=IN">IN</a>);
    <span class='reserved'>my</span> @<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = <span class='reserved'>split</span>;
    <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>[0];
}

<span class="comment"># Determine if a given pre-formatted BLAST database file is part of a</span>
<span class="comment"># multi-volume database</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=is_multivolume_db">is_multivolume_db</a>
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>return</span> 1 <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> =~ /\.\<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=d">d</a>{2,3}\.<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=tar">tar</a>\.gz$/);
    <span class='reserved'>return</span> 0;
}

<span class="comment"># Extracts the database name from the pre-formatted BLAST database archive file</span>
<span class="comment"># name</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=extract_db_name">extract_db_name</a>
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = <span class="string">""</span>;
    <span class='reserved'>if</span> (&amp;<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=is_multivolume_db">is_multivolume_db</a>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a>)) {
        $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = $1 <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> =~ <span class="extra">m/(.*)\.\d{2,3}\.tar\.gz$/</span>);
    } <span class='reserved'>else</span> {
        $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = $1 <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=file">file</a> =~ <span class="extra">m/(.*)\.tar\.gz$/</span>);
    }
    <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>;
}

<span class="comment"># Returns the number of volumes for a BLAST database given the file name of a</span>
<span class="comment"># pre-formatted BLAST database and the list of all databases to download</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_num_volumes">get_num_volumes</a>
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = 0;
    <span class='reserved'>foreach</span> (@<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=_">_</a>) {
        <span class='reserved'>if</span> (/$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=db">db</a>/) {
            <span class='reserved'>if</span> (/.*\.(\<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=d">d</a>{2,3})\.<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=tar">tar</a>\.gz$/) {
                $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = <span class='reserved'>int</span>($1) <span class='reserved'>if</span> (<span class='reserved'>int</span>($1) &gt; $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>);
            }
        }
    }
    <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> + 1;
}

<span class="comment"># Retrieves the name of the 'subdirectory' where the latest BLASTDBs residue in GCP</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_latest_dir">get_latest_dir</a>
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=source">source</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a> = <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCS_URL">GCS_URL</a> . <span class="string">"/"</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCP_BUCKET">GCP_BUCKET</a> . <span class="string">"/latest-dir"</span>;
    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a> = <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AWS_URL">AWS_URL</a> . <span class="string">"/"</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AWS_BUCKET">AWS_BUCKET</a> . <span class="string">"/latest-dir"</span> <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=source">source</a> <span class='reserved'>eq</span> <span class="string">"AWS"</span>);
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> = <span class="string">"$curl -s $url"</span>;
    <span class='reserved'>print</span> <span class="string">"$cmd\n"</span> <span class='reserved'>if</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>;
    <span class='reserved'>chomp</span>(<span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = <span class="stringx">`$cmd`</span>);
    <span class='reserved'>unless</span> (<span class='reserved'>length</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>)) {
        <span class='reserved'>print</span> STDERR <span class="string">"ERROR: Missing file $url, please try again or report to blast-help\@ncbi.nlm.nih.gov\n"</span>;
        <span class='reserved'>exit</span>(<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=EXIT_FAILURE">EXIT_FAILURE</a>);
    }
    <span class='reserved'>print</span> <span class="string">"$source latest-dir: '$retval'\n"</span> <span class='reserved'>if</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>;
    <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>;
}

<span class="comment"># Fetches the JSON text containing the BLASTDB metadata in GCS</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_blastdb_metadata">get_blastdb_metadata</a>
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=source">source</a> = <span class='reserved'>shift</span>;
    <span class='reserved'>my</span> $latest_dir = <span class='reserved'>shift</span>;
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a> = <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCS_URL">GCS_URL</a> . <span class="string">"/"</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=GCP_BUCKET">GCP_BUCKET</a> . <span class="string">"/$latest_dir/"</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=BLASTDB_MANIFEST">BLASTDB_MANIFEST</a>;
    $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a> = <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AWS_URL">AWS_URL</a> . <span class="string">"/"</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=AWS_BUCKET">AWS_BUCKET</a> . <span class="string">"/$latest_dir/"</span> . <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=BLASTDB_MANIFEST">BLASTDB_MANIFEST</a> <span class='reserved'>if</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=source">source</a> <span class='reserved'>eq</span> <span class="string">"AWS"</span>);
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=cmd">cmd</a> = <span class="string">"curl -sf $url"</span>;
    <span class='reserved'>print</span> <span class="string">"$cmd\n"</span> <span class='reserved'>if</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=DEBUG">DEBUG</a>;
    <span class='reserved'>chomp</span>(<span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = <span class="stringx">`$cmd`</span>);
    <span class='reserved'>return</span> ($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>, $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=url">url</a>);
}

<span class="comment"># Returns the path to the gsutil utility or undef if it is not found</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_gsutil_path">get_gsutil_path</a>
{
    <span class='reserved'>foreach</span> (<span class="extra">qw(/google/google-cloud-sdk/bin /usr/local/bin /usr/bin /snap/bin)</span>) {
        <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a> = <span class="string">"$_/gsutil"</span>;
        <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a> <span class='reserved'>if</span> (<span class='reserved'>-f</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a>);
    }
    <span class='reserved'>return</span> <span class='reserved'>undef</span>;
}

<span class="comment"># Returns the path to the aws CLI utility or undef if it is not found</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_awscli_path">get_awscli_path</a>
{
    <span class='reserved'>foreach</span> (<span class="extra">qw(/usr/local/bin /usr/bin)</span>) {
        <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a> = <span class="string">"$_/aws"</span>;
        <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a> <span class='reserved'>if</span> (<span class='reserved'>-f</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a>);
    }
    <span class='reserved'>return</span> <span class='reserved'>undef</span>;
}

<span class="comment"># Returns the number of cores, or 1 if unknown</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_num_cores">get_num_cores</a>
{
    <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = 1;
    <span class='reserved'>if</span> ($^O =~ /linux/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
        <span class='reserved'>open</span> <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=fh">fh</a>, <span class="string">"/proc/cpuinfo"</span> <span class='reserved'>or</span> <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>;
        $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = <span class='reserved'>scalar</span>(<span class='reserved'>map</span> /^<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=processor">processor</a>/, &lt;$<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=fh">fh</a>&gt;);
        <span class='reserved'>close</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=fh">fh</a>);
    } <span class='reserved'>elsif</span> ($^O =~ /darwin/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
        <span class='reserved'>chomp</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = <span class="stringx">`/usr/sbin/sysctl -n hw.ncpu`</span>);
    }
    <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>;
}

<span class="comment"># Returns the path to the curl utility, or undef if it is not found</span>
<span class='reserved'>sub</span> <a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=get_curl_path">get_curl_path</a>
{
    <span class='reserved'>foreach</span> (<span class="extra">qw(/usr/local/bin /usr/bin)</span>) {
        <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a> = <span class="string">"$_/curl"</span>;
        <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a> <span class='reserved'>if</span> (<span class='reserved'>-f</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=path">path</a>);
    }
    <span class='reserved'>if</span> ($^O =~ /mswin/<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=i">i</a>) {
        <span class='reserved'>my</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> = <span class="stringx">`where curl`</span>;
        <span class='reserved'>if</span> (<span class='reserved'>defined</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>) {
            <span class='reserved'>chomp</span>($<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>);
            <span class='reserved'>return</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a> <span class='reserved'>if</span> (<span class='reserved'>-f</span> $<a class='fid' href="/IEB/ToolBox/CPP_DOC/lxr/ident?i=retval">retval</a>);
        }
    }
    <span class='reserved'>return</span> <span class='reserved'>undef</span>;
}

<span class='reserved'>__END__</span>

<span class="comment-d">=head1 NAME</span>
<span class="comment-d"></span>
<span class="comment-d">B&lt;update_blastdb.pl&gt; - Download pre-formatted BLAST databases</span>
<span class="comment-d"></span>
<span class="comment-d">=head1 SYNOPSIS</span>
<span class="comment-d"></span>
<span class="comment-d">update_blastdb.pl [options] blastdb ...</span>
<span class="comment-d"></span>
<span class="comment-d">=head1 OPTIONS</span>
<span class="comment-d"></span>
<span class="comment-d">=over 2</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--source&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Location to download BLAST databases from (default: auto-detect closest location).</span>
<span class="comment-d">Supported values: ncbi, aws, or gcp.</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--decompress&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Downloads, decompresses the archives in the current working directory, and</span>
<span class="comment-d">deletes the downloaded archive to save disk space, while preserving the</span>
<span class="comment-d">archive checksum files (default: false).</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--showall&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Show all available pre-formatted BLAST databases (default: false). The output</span>
<span class="comment-d">of this option lists the database names which should be used when</span>
<span class="comment-d">requesting downloads or updates using this script.</span>
<span class="comment-d"></span>
<span class="comment-d">It accepts the optional arguments: 'tsv' and 'pretty' to produce tab-separated values</span>
<span class="comment-d">and a human-readable format respectively. These parameters elicit the display of</span>
<span class="comment-d">additional metadata if this is available to the program.</span>
<span class="comment-d">This metadata is displayed in columnar format; the columns represent:</span>
<span class="comment-d"></span>
<span class="comment-d">name, description, size in gigabytes, date of last update (YYYY-MM-DD format).</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--blastdb_version&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Specify which BLAST database version to download (default: 5).</span>
<span class="comment-d">Supported values: 4, 5</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--passive&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Use passive FTP, useful when behind a firewall or working in the cloud (default: true).</span>
<span class="comment-d">To disable passive FTP, configure this option as follows: --passive no</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--timeout&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Timeout on connection to NCBI (default: 120 seconds).</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--force&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Force download even if there is a archive already on local directory (default:</span>
<span class="comment-d">false).</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--verbose&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Increment verbosity level (default: 1). Repeat this option multiple times to </span>
<span class="comment-d">increase the verbosity level (maximum 2).</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--quiet&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Produce no output (default: false). Overrides the B&lt;--verbose&gt; option.</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--version&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Prints this script's version. Overrides all other options.</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--num_threads&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Sets the number of threads to utilize to perform downloads in parallel when data comes from the cloud.</span>
<span class="comment-d">Defaults to use all available CPUs on the machine (Linux and macos only).</span>
<span class="comment-d"></span>
<span class="comment-d">=item B&lt;--legacy_exit_code&gt;</span>
<span class="comment-d"></span>
<span class="comment-d">Enables exit codes from prior to version 581818, BLAST+ 2.10.0 release, for</span>
<span class="comment-d">downloads from NCBI only. This option is meant to be used by legacy applications that rely</span>
<span class="comment-d">on this exit codes:</span>
<span class="comment-d">0 for successful operations that result in no downloads, 1 for successful</span>
<span class="comment-d">downloads, and 2 for errors.</span>
<span class="comment-d"></span>
<span class="comment-d">=back</span>
<span class="comment-d"></span>
<span class="comment-d">=head1 DESCRIPTION</span>
<span class="comment-d"></span>
<span class="comment-d">This script will download the pre-formatted BLAST databases requested in the</span>
<span class="comment-d">command line from the NCBI ftp site.</span>
<span class="comment-d"></span>
<span class="comment-d">=head1 EXIT CODES</span>
<span class="comment-d"></span>
<span class="comment-d">This script returns 0 on successful operations and non-zero on errors.</span>
<span class="comment-d"></span>
<span class="comment-d">=head1 DEPENDENCIES</span>
<span class="comment-d"></span>
<span class="comment-d">This script depends on curl for retrieval from cloud service providers.</span>
<span class="comment-d"></span>
<span class="comment-d">=head1 BUGS</span>
<span class="comment-d"></span>
<span class="comment-d">Please report them to <a class='offshore' href="mailto:blast-help@ncbi.nlm.nih.gov">&lt;blast-help@ncbi.nlm.nih.gov&gt;</a></span>
<span class="comment-d"></span>
<span class="comment-d">=head1 COPYRIGHT</span>
<span class="comment-d"></span>
<span class="comment-d">See PUBLIC DOMAIN NOTICE included at the top of this script.</span>
<span class="comment-d"></span>
<span class="comment-d">=cut</span></pre>
<pre class="filecontent-num">
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0001" name="0001">0001</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0002" name="0002">0002</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0003" name="0003">0003</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0004" name="0004">0004</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0005" name="0005">0005</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0006" name="0006">0006</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0007" name="0007">0007</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0008" name="0008">0008</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0009" name="0009">0009</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0010" name="0010">0010</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0011" name="0011">0011</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0012" name="0012">0012</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0013" name="0013">0013</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0014" name="0014">0014</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0015" name="0015">0015</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0016" name="0016">0016</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0017" name="0017">0017</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0018" name="0018">0018</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0019" name="0019">0019</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0020" name="0020">0020</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0021" name="0021">0021</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0022" name="0022">0022</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0023" name="0023">0023</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0024" name="0024">0024</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0025" name="0025">0025</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0026" name="0026">0026</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0027" name="0027">0027</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0028" name="0028">0028</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0029" name="0029">0029</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0030" name="0030">0030</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0031" name="0031">0031</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0032" name="0032">0032</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0033" name="0033">0033</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0034" name="0034">0034</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0035" name="0035">0035</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0036" name="0036">0036</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0037" name="0037">0037</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0038" name="0038">0038</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0039" name="0039">0039</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0040" name="0040">0040</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0041" name="0041">0041</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0042" name="0042">0042</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0043" name="0043">0043</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0044" name="0044">0044</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0045" name="0045">0045</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0046" name="0046">0046</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0047" name="0047">0047</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0048" name="0048">0048</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0049" name="0049">0049</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0050" name="0050">0050</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0051" name="0051">0051</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0052" name="0052">0052</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0053" name="0053">0053</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0054" name="0054">0054</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0055" name="0055">0055</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0056" name="0056">0056</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0057" name="0057">0057</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0058" name="0058">0058</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0059" name="0059">0059</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0060" name="0060">0060</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0061" name="0061">0061</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0062" name="0062">0062</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0063" name="0063">0063</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0064" name="0064">0064</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0065" name="0065">0065</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0066" name="0066">0066</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0067" name="0067">0067</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0068" name="0068">0068</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0069" name="0069">0069</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0070" name="0070">0070</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0071" name="0071">0071</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0072" name="0072">0072</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0073" name="0073">0073</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0074" name="0074">0074</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0075" name="0075">0075</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0076" name="0076">0076</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0077" name="0077">0077</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0078" name="0078">0078</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0079" name="0079">0079</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0080" name="0080">0080</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0081" name="0081">0081</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0082" name="0082">0082</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0083" name="0083">0083</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0084" name="0084">0084</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0085" name="0085">0085</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0086" name="0086">0086</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0087" name="0087">0087</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0088" name="0088">0088</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0089" name="0089">0089</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0090" name="0090">0090</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0091" name="0091">0091</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0092" name="0092">0092</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0093" name="0093">0093</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0094" name="0094">0094</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0095" name="0095">0095</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0096" name="0096">0096</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0097" name="0097">0097</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0098" name="0098">0098</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0099" name="0099">0099</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0100" name="0100">0100</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0101" name="0101">0101</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0102" name="0102">0102</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0103" name="0103">0103</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0104" name="0104">0104</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0105" name="0105">0105</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0106" name="0106">0106</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0107" name="0107">0107</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0108" name="0108">0108</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0109" name="0109">0109</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0110" name="0110">0110</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0111" name="0111">0111</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0112" name="0112">0112</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0113" name="0113">0113</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0114" name="0114">0114</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0115" name="0115">0115</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0116" name="0116">0116</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0117" name="0117">0117</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0118" name="0118">0118</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0119" name="0119">0119</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0120" name="0120">0120</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0121" name="0121">0121</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0122" name="0122">0122</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0123" name="0123">0123</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0124" name="0124">0124</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0125" name="0125">0125</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0126" name="0126">0126</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0127" name="0127">0127</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0128" name="0128">0128</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0129" name="0129">0129</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0130" name="0130">0130</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0131" name="0131">0131</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0132" name="0132">0132</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0133" name="0133">0133</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0134" name="0134">0134</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0135" name="0135">0135</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0136" name="0136">0136</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0137" name="0137">0137</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0138" name="0138">0138</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0139" name="0139">0139</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0140" name="0140">0140</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0141" name="0141">0141</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0142" name="0142">0142</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0143" name="0143">0143</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0144" name="0144">0144</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0145" name="0145">0145</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0146" name="0146">0146</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0147" name="0147">0147</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0148" name="0148">0148</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0149" name="0149">0149</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0150" name="0150">0150</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0151" name="0151">0151</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0152" name="0152">0152</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0153" name="0153">0153</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0154" name="0154">0154</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0155" name="0155">0155</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0156" name="0156">0156</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0157" name="0157">0157</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0158" name="0158">0158</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0159" name="0159">0159</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0160" name="0160">0160</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0161" name="0161">0161</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0162" name="0162">0162</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0163" name="0163">0163</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0164" name="0164">0164</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0165" name="0165">0165</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0166" name="0166">0166</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0167" name="0167">0167</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0168" name="0168">0168</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0169" name="0169">0169</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0170" name="0170">0170</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0171" name="0171">0171</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0172" name="0172">0172</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0173" name="0173">0173</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0174" name="0174">0174</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0175" name="0175">0175</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0176" name="0176">0176</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0177" name="0177">0177</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0178" name="0178">0178</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0179" name="0179">0179</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0180" name="0180">0180</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0181" name="0181">0181</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0182" name="0182">0182</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0183" name="0183">0183</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0184" name="0184">0184</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0185" name="0185">0185</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0186" name="0186">0186</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0187" name="0187">0187</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0188" name="0188">0188</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0189" name="0189">0189</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0190" name="0190">0190</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0191" name="0191">0191</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0192" name="0192">0192</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0193" name="0193">0193</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0194" name="0194">0194</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0195" name="0195">0195</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0196" name="0196">0196</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0197" name="0197">0197</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0198" name="0198">0198</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0199" name="0199">0199</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0200" name="0200">0200</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0201" name="0201">0201</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0202" name="0202">0202</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0203" name="0203">0203</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0204" name="0204">0204</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0205" name="0205">0205</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0206" name="0206">0206</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0207" name="0207">0207</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0208" name="0208">0208</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0209" name="0209">0209</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0210" name="0210">0210</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0211" name="0211">0211</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0212" name="0212">0212</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0213" name="0213">0213</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0214" name="0214">0214</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0215" name="0215">0215</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0216" name="0216">0216</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0217" name="0217">0217</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0218" name="0218">0218</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0219" name="0219">0219</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0220" name="0220">0220</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0221" name="0221">0221</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0222" name="0222">0222</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0223" name="0223">0223</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0224" name="0224">0224</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0225" name="0225">0225</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0226" name="0226">0226</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0227" name="0227">0227</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0228" name="0228">0228</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0229" name="0229">0229</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0230" name="0230">0230</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0231" name="0231">0231</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0232" name="0232">0232</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0233" name="0233">0233</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0234" name="0234">0234</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0235" name="0235">0235</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0236" name="0236">0236</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0237" name="0237">0237</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0238" name="0238">0238</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0239" name="0239">0239</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0240" name="0240">0240</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0241" name="0241">0241</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0242" name="0242">0242</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0243" name="0243">0243</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0244" name="0244">0244</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0245" name="0245">0245</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0246" name="0246">0246</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0247" name="0247">0247</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0248" name="0248">0248</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0249" name="0249">0249</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0250" name="0250">0250</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0251" name="0251">0251</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0252" name="0252">0252</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0253" name="0253">0253</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0254" name="0254">0254</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0255" name="0255">0255</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0256" name="0256">0256</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0257" name="0257">0257</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0258" name="0258">0258</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0259" name="0259">0259</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0260" name="0260">0260</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0261" name="0261">0261</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0262" name="0262">0262</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0263" name="0263">0263</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0264" name="0264">0264</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0265" name="0265">0265</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0266" name="0266">0266</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0267" name="0267">0267</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0268" name="0268">0268</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0269" name="0269">0269</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0270" name="0270">0270</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0271" name="0271">0271</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0272" name="0272">0272</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0273" name="0273">0273</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0274" name="0274">0274</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0275" name="0275">0275</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0276" name="0276">0276</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0277" name="0277">0277</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0278" name="0278">0278</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0279" name="0279">0279</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0280" name="0280">0280</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0281" name="0281">0281</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0282" name="0282">0282</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0283" name="0283">0283</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0284" name="0284">0284</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0285" name="0285">0285</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0286" name="0286">0286</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0287" name="0287">0287</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0288" name="0288">0288</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0289" name="0289">0289</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0290" name="0290">0290</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0291" name="0291">0291</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0292" name="0292">0292</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0293" name="0293">0293</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0294" name="0294">0294</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0295" name="0295">0295</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0296" name="0296">0296</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0297" name="0297">0297</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0298" name="0298">0298</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0299" name="0299">0299</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0300" name="0300">0300</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0301" name="0301">0301</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0302" name="0302">0302</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0303" name="0303">0303</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0304" name="0304">0304</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0305" name="0305">0305</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0306" name="0306">0306</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0307" name="0307">0307</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0308" name="0308">0308</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0309" name="0309">0309</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0310" name="0310">0310</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0311" name="0311">0311</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0312" name="0312">0312</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0313" name="0313">0313</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0314" name="0314">0314</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0315" name="0315">0315</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0316" name="0316">0316</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0317" name="0317">0317</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0318" name="0318">0318</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0319" name="0319">0319</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0320" name="0320">0320</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0321" name="0321">0321</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0322" name="0322">0322</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0323" name="0323">0323</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0324" name="0324">0324</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0325" name="0325">0325</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0326" name="0326">0326</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0327" name="0327">0327</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0328" name="0328">0328</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0329" name="0329">0329</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0330" name="0330">0330</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0331" name="0331">0331</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0332" name="0332">0332</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0333" name="0333">0333</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0334" name="0334">0334</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0335" name="0335">0335</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0336" name="0336">0336</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0337" name="0337">0337</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0338" name="0338">0338</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0339" name="0339">0339</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0340" name="0340">0340</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0341" name="0341">0341</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0342" name="0342">0342</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0343" name="0343">0343</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0344" name="0344">0344</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0345" name="0345">0345</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0346" name="0346">0346</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0347" name="0347">0347</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0348" name="0348">0348</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0349" name="0349">0349</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0350" name="0350">0350</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0351" name="0351">0351</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0352" name="0352">0352</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0353" name="0353">0353</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0354" name="0354">0354</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0355" name="0355">0355</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0356" name="0356">0356</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0357" name="0357">0357</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0358" name="0358">0358</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0359" name="0359">0359</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0360" name="0360">0360</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0361" name="0361">0361</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0362" name="0362">0362</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0363" name="0363">0363</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0364" name="0364">0364</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0365" name="0365">0365</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0366" name="0366">0366</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0367" name="0367">0367</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0368" name="0368">0368</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0369" name="0369">0369</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0370" name="0370">0370</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0371" name="0371">0371</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0372" name="0372">0372</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0373" name="0373">0373</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0374" name="0374">0374</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0375" name="0375">0375</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0376" name="0376">0376</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0377" name="0377">0377</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0378" name="0378">0378</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0379" name="0379">0379</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0380" name="0380">0380</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0381" name="0381">0381</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0382" name="0382">0382</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0383" name="0383">0383</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0384" name="0384">0384</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0385" name="0385">0385</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0386" name="0386">0386</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0387" name="0387">0387</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0388" name="0388">0388</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0389" name="0389">0389</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0390" name="0390">0390</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0391" name="0391">0391</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0392" name="0392">0392</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0393" name="0393">0393</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0394" name="0394">0394</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0395" name="0395">0395</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0396" name="0396">0396</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0397" name="0397">0397</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0398" name="0398">0398</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0399" name="0399">0399</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0400" name="0400">0400</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0401" name="0401">0401</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0402" name="0402">0402</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0403" name="0403">0403</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0404" name="0404">0404</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0405" name="0405">0405</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0406" name="0406">0406</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0407" name="0407">0407</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0408" name="0408">0408</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0409" name="0409">0409</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0410" name="0410">0410</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0411" name="0411">0411</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0412" name="0412">0412</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0413" name="0413">0413</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0414" name="0414">0414</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0415" name="0415">0415</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0416" name="0416">0416</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0417" name="0417">0417</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0418" name="0418">0418</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0419" name="0419">0419</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0420" name="0420">0420</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0421" name="0421">0421</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0422" name="0422">0422</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0423" name="0423">0423</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0424" name="0424">0424</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0425" name="0425">0425</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0426" name="0426">0426</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0427" name="0427">0427</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0428" name="0428">0428</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0429" name="0429">0429</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0430" name="0430">0430</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0431" name="0431">0431</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0432" name="0432">0432</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0433" name="0433">0433</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0434" name="0434">0434</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0435" name="0435">0435</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0436" name="0436">0436</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0437" name="0437">0437</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0438" name="0438">0438</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0439" name="0439">0439</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0440" name="0440">0440</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0441" name="0441">0441</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0442" name="0442">0442</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0443" name="0443">0443</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0444" name="0444">0444</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0445" name="0445">0445</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0446" name="0446">0446</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0447" name="0447">0447</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0448" name="0448">0448</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0449" name="0449">0449</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0450" name="0450">0450</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0451" name="0451">0451</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0452" name="0452">0452</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0453" name="0453">0453</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0454" name="0454">0454</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0455" name="0455">0455</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0456" name="0456">0456</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0457" name="0457">0457</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0458" name="0458">0458</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0459" name="0459">0459</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0460" name="0460">0460</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0461" name="0461">0461</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0462" name="0462">0462</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0463" name="0463">0463</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0464" name="0464">0464</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0465" name="0465">0465</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0466" name="0466">0466</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0467" name="0467">0467</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0468" name="0468">0468</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0469" name="0469">0469</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0470" name="0470">0470</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0471" name="0471">0471</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0472" name="0472">0472</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0473" name="0473">0473</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0474" name="0474">0474</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0475" name="0475">0475</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0476" name="0476">0476</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0477" name="0477">0477</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0478" name="0478">0478</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0479" name="0479">0479</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0480" name="0480">0480</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0481" name="0481">0481</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0482" name="0482">0482</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0483" name="0483">0483</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0484" name="0484">0484</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0485" name="0485">0485</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0486" name="0486">0486</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0487" name="0487">0487</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0488" name="0488">0488</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0489" name="0489">0489</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0490" name="0490">0490</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0491" name="0491">0491</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0492" name="0492">0492</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0493" name="0493">0493</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0494" name="0494">0494</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0495" name="0495">0495</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0496" name="0496">0496</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0497" name="0497">0497</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0498" name="0498">0498</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0499" name="0499">0499</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0500" name="0500">0500</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0501" name="0501">0501</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0502" name="0502">0502</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0503" name="0503">0503</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0504" name="0504">0504</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0505" name="0505">0505</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0506" name="0506">0506</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0507" name="0507">0507</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0508" name="0508">0508</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0509" name="0509">0509</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0510" name="0510">0510</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0511" name="0511">0511</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0512" name="0512">0512</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0513" name="0513">0513</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0514" name="0514">0514</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0515" name="0515">0515</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0516" name="0516">0516</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0517" name="0517">0517</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0518" name="0518">0518</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0519" name="0519">0519</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0520" name="0520">0520</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0521" name="0521">0521</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0522" name="0522">0522</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0523" name="0523">0523</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0524" name="0524">0524</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0525" name="0525">0525</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0526" name="0526">0526</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0527" name="0527">0527</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0528" name="0528">0528</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0529" name="0529">0529</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0530" name="0530">0530</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0531" name="0531">0531</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0532" name="0532">0532</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0533" name="0533">0533</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0534" name="0534">0534</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0535" name="0535">0535</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0536" name="0536">0536</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0537" name="0537">0537</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0538" name="0538">0538</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0539" name="0539">0539</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0540" name="0540">0540</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0541" name="0541">0541</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0542" name="0542">0542</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0543" name="0543">0543</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0544" name="0544">0544</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0545" name="0545">0545</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0546" name="0546">0546</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0547" name="0547">0547</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0548" name="0548">0548</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0549" name="0549">0549</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0550" name="0550">0550</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0551" name="0551">0551</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0552" name="0552">0552</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0553" name="0553">0553</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0554" name="0554">0554</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0555" name="0555">0555</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0556" name="0556">0556</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0557" name="0557">0557</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0558" name="0558">0558</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0559" name="0559">0559</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0560" name="0560">0560</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0561" name="0561">0561</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0562" name="0562">0562</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0563" name="0563">0563</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0564" name="0564">0564</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0565" name="0565">0565</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0566" name="0566">0566</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0567" name="0567">0567</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0568" name="0568">0568</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0569" name="0569">0569</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0570" name="0570">0570</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0571" name="0571">0571</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0572" name="0572">0572</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0573" name="0573">0573</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0574" name="0574">0574</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0575" name="0575">0575</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0576" name="0576">0576</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0577" name="0577">0577</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0578" name="0578">0578</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0579" name="0579">0579</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0580" name="0580">0580</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0581" name="0581">0581</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0582" name="0582">0582</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0583" name="0583">0583</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0584" name="0584">0584</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0585" name="0585">0585</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0586" name="0586">0586</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0587" name="0587">0587</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0588" name="0588">0588</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0589" name="0589">0589</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0590" name="0590">0590</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0591" name="0591">0591</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0592" name="0592">0592</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0593" name="0593">0593</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0594" name="0594">0594</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0595" name="0595">0595</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0596" name="0596">0596</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0597" name="0597">0597</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0598" name="0598">0598</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0599" name="0599">0599</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0600" name="0600">0600</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0601" name="0601">0601</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0602" name="0602">0602</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0603" name="0603">0603</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0604" name="0604">0604</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0605" name="0605">0605</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0606" name="0606">0606</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0607" name="0607">0607</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0608" name="0608">0608</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0609" name="0609">0609</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0610" name="0610">0610</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0611" name="0611">0611</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0612" name="0612">0612</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0613" name="0613">0613</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0614" name="0614">0614</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0615" name="0615">0615</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0616" name="0616">0616</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0617" name="0617">0617</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0618" name="0618">0618</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0619" name="0619">0619</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0620" name="0620">0620</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0621" name="0621">0621</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0622" name="0622">0622</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0623" name="0623">0623</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0624" name="0624">0624</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0625" name="0625">0625</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0626" name="0626">0626</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0627" name="0627">0627</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0628" name="0628">0628</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0629" name="0629">0629</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0630" name="0630">0630</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0631" name="0631">0631</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0632" name="0632">0632</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0633" name="0633">0633</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0634" name="0634">0634</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0635" name="0635">0635</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0636" name="0636">0636</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0637" name="0637">0637</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0638" name="0638">0638</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0639" name="0639">0639</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0640" name="0640">0640</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0641" name="0641">0641</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0642" name="0642">0642</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0643" name="0643">0643</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0644" name="0644">0644</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0645" name="0645">0645</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0646" name="0646">0646</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0647" name="0647">0647</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0648" name="0648">0648</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0649" name="0649">0649</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0650" name="0650">0650</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0651" name="0651">0651</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0652" name="0652">0652</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0653" name="0653">0653</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0654" name="0654">0654</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0655" name="0655">0655</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0656" name="0656">0656</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0657" name="0657">0657</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0658" name="0658">0658</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0659" name="0659">0659</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0660" name="0660">0660</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0661" name="0661">0661</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0662" name="0662">0662</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0663" name="0663">0663</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0664" name="0664">0664</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0665" name="0665">0665</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0666" name="0666">0666</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0667" name="0667">0667</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0668" name="0668">0668</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0669" name="0669">0669</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0670" name="0670">0670</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0671" name="0671">0671</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0672" name="0672">0672</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0673" name="0673">0673</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0674" name="0674">0674</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0675" name="0675">0675</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0676" name="0676">0676</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0677" name="0677">0677</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0678" name="0678">0678</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0679" name="0679">0679</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0680" name="0680">0680</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0681" name="0681">0681</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0682" name="0682">0682</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0683" name="0683">0683</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0684" name="0684">0684</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0685" name="0685">0685</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0686" name="0686">0686</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0687" name="0687">0687</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0688" name="0688">0688</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0689" name="0689">0689</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0690" name="0690">0690</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0691" name="0691">0691</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0692" name="0692">0692</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0693" name="0693">0693</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0694" name="0694">0694</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0695" name="0695">0695</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0696" name="0696">0696</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0697" name="0697">0697</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0698" name="0698">0698</a>
<a class="fline" href="/IEB/ToolBox/CPP_DOC/lxr/source/src/app/blast/update_blastdb.pl#0699" name="0699">0699</a>
</pre>
</div>
</main>
<footer>
<hr>
<div align="center">
     [&nbsp;<span class='modes-sel'>Source navigation</span>&nbsp;]&nbsp;&nbsp;  [&nbsp;<a class='modes' href="/IEB/ToolBox/CPP_DOC/lxr/diff/src/app/blast/update_blastdb.pl">Diff markup</a>&nbsp;]&nbsp;&nbsp;  [&nbsp;<a class='modes' href="/IEB/ToolBox/CPP_DOC/lxr/ident">Identifier search</a>&nbsp;]&nbsp;&nbsp;  [&nbsp;<a class="modes" href="/IEB/ToolBox/CPP_DOC/lxr/search">General search</a>&nbsp;]&nbsp;&nbsp; 
</div>
<hr>
<div class="version">
    Generated by the LXR 2.3.5.&nbsp;&mdash;&nbsp;<span class=indexstate>Indexed on 2021-01-24 23:38:05 UTC</span>

</div>
</footer>
</body>
</html>
