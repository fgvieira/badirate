# NOTE: Derived from blib/lib/Inline.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Inline;

#line 1006 "blib/lib/Inline.pm (autosplit into blib/lib/auto/Inline/env_untaint.al)"
#==============================================================================
# Blindly untaint tainted fields in %ENV.
#==============================================================================
sub env_untaint {
    my $o = shift;
    warn "In Inline::env_untaint() : Blindly untainting tainted fields in %ENV.\n" unless $o->{CONFIG}{NO_UNTAINT_WARN};

    for (keys %ENV) {
	($ENV{$_}) = $ENV{$_} =~ /(.*)/;
    }

    $ENV{PATH} = $^O eq 'MSWin32' ?
                 join ';', grep {not /^\./ and -d $_
				  } split /;/, $ENV{PATH}
                 :
                 join ':', grep {not /^\./ and -d $_ and
				      not ((stat($_))[2] & 0022)
				  } split /:/, $ENV{PATH};
    map {($_) = /(.*)/} @INC;
}

# end of Inline::env_untaint
1;
