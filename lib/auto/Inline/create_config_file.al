# NOTE: Derived from blib/lib/Inline.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Inline;

#line 749 "blib/lib/Inline.pm (autosplit into blib/lib/auto/Inline/create_config_file.al)"
#==============================================================================
# Auto-detect installed Inline language support modules
#==============================================================================
sub create_config_file {
    my ($o, $dir) = @_;

    # This subroutine actually fires off another instance of perl.
    # with arguments that make this routine get called again.
    # That way the queried modules don't stay loaded.
    if (defined $o) {
	($dir) = $dir =~ /(.*)/s if UNTAINT;
	my $perl = $Config{perlpath};
        $perl = $^X unless -f $perl;
	($perl) = $perl =~ /(.*)/s if UNTAINT;
	local $ENV{PERL5LIB} if defined $ENV{PERL5LIB};
	local $ENV{PERL5OPT} if defined $ENV{PERL5OPT};
	my $inline = $INC{'Inline.pm'};
        $inline ||= File::Spec->curdir();
        my($v,$d,$f) = File::Spec->splitpath($inline);
        $f = "" if $f eq 'Inline.pm';
        $inline = File::Spec->catpath($v,$d,$f);
        #my $INC = "-I$inline -I" .
        #          join(" -I", grep {(-d File::Spec->catdir($_,"Inline") or
        #                             -d File::Spec->catdir($_,"auto","Inline")
	  #		            )} @INC);
	#system "$perl $INC -MInline=_CONFIG_ -e1 $dir"
        my @INC = map { "-I$_" }
       ($inline,
        grep {(-d File::Spec->catdir($_,"Inline") or -d File::Spec->catdir($_,"auto","Inline"))} @INC);
       system $perl, @INC, "-MInline=_CONFIG_", "-e1", "$dir"
	  and croak M20_config_creation_failed($dir);
	return;
    }

    my ($lib, $mod, $register, %checked,
	%languages, %types, %modules, %suffixes);
  LIB:
    for my $lib (@INC) {
        next unless -d File::Spec->catdir($lib,"Inline");
        opendir LIB, File::Spec->catdir($lib,"Inline")
          or warn(M21_opendir_failed(File::Spec->catdir($lib,"Inline"))), next;
	while ($mod = readdir(LIB)) {
	    next unless $mod =~ /\.pm$/;
	    $mod =~ s/\.pm$//;
	    next LIB if ($checked{$mod}++);
	    if ($mod eq 'Config') {     # Skip Inline::Config
		warn M14_usage_Config();
		next;
	    }
	    next if $mod =~ /^(MakeMaker|denter|messages)$/;
	    eval "require Inline::$mod;";
            warn($@), next if $@;
	    eval "\$register=&Inline::${mod}::register";
	    next if $@;
	    my $language = ($register->{language})
	      or warn(M22_usage_register($mod)), next;
	    for (@{$register->{aliases}}) {
		warn(M23_usage_alias_used($mod, $_, $languages{$_})), next
		  if defined $languages{$_};
		$languages{$_} = $language;
	    }
	    $languages{$language} = $language;
	    $types{$language} = $register->{type};
	    $modules{$language} = "Inline::$mod";
	    $suffixes{$language} = $register->{suffix};
	}
	closedir LIB;
    }

    my $file = File::Spec->catfile($ARGV[0],"config");
    open CONFIG, "> $file" or croak M24_open_for_output_failed($file);
    print CONFIG Inline::denter->new()
      ->indent(*version => $Inline::VERSION,
	       *languages => \%languages,
	       *types => \%types,
	       *modules => \%modules,
	       *suffixes => \%suffixes,
	      );
    close CONFIG;
    exit 0;
}

# end of Inline::create_config_file
1;
