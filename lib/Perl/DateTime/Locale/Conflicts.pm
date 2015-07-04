package # hide from PAUSE
    DateTime::Locale::Conflicts;

use strict;
use warnings;

# this module was generated with Dist::Zilla::Plugin::Conflicts 0.17

use Dist::CheckConflicts
    -dist      => 'DateTime::Locale',
    -conflicts => {
        'DateTime::Format::Strptime' => '1.1000',
    },
    -also => [ qw(
        Carp
        Dist::CheckConflicts
        List::MoreUtils
        Params::Validate
        base
        strict
        utf8
        warnings
    ) ],

;

1;

# ABSTRACT: Provide information on conflicts for DateTime::Locale
# Dist::Zilla: -PodWeaver
