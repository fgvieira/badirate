# NOTE: Derived from blib/lib/Inline.pm.
# Changes made here will be lost when autosplit is run again.
# See AutoSplit.pm.
package Inline;

#line 1539 "blib/lib/Inline.pm (autosplit into blib/lib/auto/Inline/M19_usage_language.al)"
sub M19_usage_language {
    my ($language, $directory) = @_;
    return <<END;
Error. You have specified '$language' as an Inline programming language.

I currently only know about the following languages:
    ${ defined $Inline::languages ?
       \ join(', ', sort keys %$Inline::languages) : \ ''
     }

If you have installed a support module for this language, try deleting the
config file from the following Inline DIRECTORY, and run again:

    $directory

END
}

# end of Inline::M19_usage_language
1;
