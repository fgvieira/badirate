# NOTE: Derived from blib/lib/Inline.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Inline;

#line 686 "blib/lib/Inline.pm (autosplit into blib/lib/auto/Inline/check_config_file.al)"
#==============================================================================
# Read the cached config file from the Inline directory. This will indicate
# whether the Language code is valid or not.
#==============================================================================
sub check_config_file {
    my ($DIRECTORY, %config);
    my $o = shift;

    croak M14_usage_Config() if %main::Inline::Config::;
    croak M63_no_source($o->{API}{pkg})
      if $o->{INLINE}{md5} eq $o->{API}{code};

    # First make sure we have the DIRECTORY
    if ($o->{CONFIG}{_INSTALL_}) {
	croak M15_usage_install_directory()
	  if $o->{CONFIG}{DIRECTORY};
	my $cwd = Cwd::cwd();
        $DIRECTORY =
          $o->{INLINE}{DIRECTORY} = File::Spec->catdir($cwd, $did);
	if (not -d $DIRECTORY) {
	    _mkdir($DIRECTORY, 0777)
	      or croak M16_DIRECTORY_mkdir_failed($DIRECTORY);
	}
    }
    else {
	$DIRECTORY = $o->{INLINE}{DIRECTORY} =
	  $o->{CONFIG}{DIRECTORY} || $o->find_temp_dir;
    }

    $o->create_config_file($DIRECTORY)
      if not -e File::Spec->catfile($DIRECTORY,"config");

    open CONFIG, "< ".File::Spec->catfile($DIRECTORY,"config")
      or croak M17_config_open_failed($DIRECTORY);
    my $config = join '', <CONFIG>;
    close CONFIG;

    croak M62_invalid_config_file(File::Spec->catfile($DIRECTORY,"config"))
      unless $config =~ /^version :/;
    ($config) = $config =~ /(.*)/s if UNTAINT;

    %config = Inline::denter->new()->undent($config);
    $Inline::languages = $config{languages};

    croak M18_error_old_version($config{version}, $DIRECTORY)
	unless (defined $config{version} and
                $config{version} =~ /TRIAL/ or
		$config{version} >= 0.40);
    croak M19_usage_language($o->{API}{language_id}, $DIRECTORY)
      unless defined $config{languages}->{$o->{API}{language_id}};
    $o->{API}{language} = $config{languages}->{$o->{API}{language_id}};
    if ($o->{API}{language} ne $o->{API}{language_id}) {
	if (defined $o->{$o->{API}{language_id}}) {
	    $o->{$o->{API}{language}} = $o->{$o->{API}{language_id}};
	    delete $o->{$o->{API}{language_id}};
	}
    }

    $o->{INLINE}{ILSM_type} = $config{types}->{$o->{API}{language}};
    $o->{INLINE}{ILSM_module} = $config{modules}->{$o->{API}{language}};
    $o->{INLINE}{ILSM_suffix} = $config{suffixes}->{$o->{API}{language}};
}

# end of Inline::check_config_file
1;
