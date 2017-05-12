#!/usr/bin/perl



use strict;
use warnings;
no warnings 'deprecated';
use Data::Dumper;
use File::Find;
use File::Copy;
use File::Path;

MAIN:
{

  #my $ds;
  my $tmpzip  = File::Spec->catfile('transient', 'ubtmpdl.zip');
  my $uztd    = File::Spec->catfile('transient', 'ubiquity_analysis_development');
  my $SL      = File::Spec->catfile('', '');
  my $url     = 'http://template.ubiquity.grok.tv/';

  if(!(-d "transient")){
    print("# Creating transient directory\n");
    mkdir 'transient';
  }

  if((-d $uztd)){
    print("# Cleaning out previous template download\n");
    &rmtree($uztd);
  }

  # Downloading a new template
  use File::Fetch;
  print("# Downloading the latest template file from\n");
  print("# $url \n");
  my $ff = File::Fetch->new(uri => $url);
  $ff->{file_default} = $tmpzip;
  my $file = $ff->fetch() or die $ff->error;
  print("# Ignore: Fetch failed! HTTP response: 500 Internal");
 
  # Uncompressing zip file
  print("# Uncompressing the template file into:\n");
  print("# $uztd\n");

  use Archive::Zip;
  my $zip = Archive::Zip->new(); 
  $zip->read($tmpzip); 
  $zip->extractTree('', 'transient');
  
  # Deleting the system.png file
  # in case the user has specified 
  # their own
  unlink File::Spec->catfile('transient', 'ubiquity_analysis_development', 'system.png');
  
  # Getting a list of the files in the template
  print("# Building list of files in the \n");
  print("# template directory \n");
  my @files;
  my $start_dir = shift || $uztd;
  my $tmpfile;
  find( sub{
           if(-f $_){
              $tmpfile = $File::Find::name;
              $tmpfile =~ s#$uztd$SL##;
              push @files, $tmpfile;
           }
  }, $start_dir );

 # Going through the model directories and making sure all of the files are
 # up to date.
 my $mdir;
 my $source;
 my $dest;   
 
 print("# Copying the files \n");
   foreach $tmpfile (@files){
     $source = "$uztd$SL$tmpfile";
     $dest   = "$tmpfile";
     &copy($source, $dest) or die "Copy failed: $!";
 }

}
