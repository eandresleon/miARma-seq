# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/8FT049ktOU/australasia.  Olson data version 2015d
#
# Do not edit this file directly.
#
package DateTime::TimeZone::Pacific::Majuro;
$DateTime::TimeZone::Pacific::Majuro::VERSION = '1.90';
use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::Pacific::Majuro::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
59958189312, #      utc_end 1900-12-31 12:35:12 (Mon)
DateTime::TimeZone::NEG_INFINITY, #  local_start
59958230400, #    local_end 1901-01-01 00:00:00 (Tue)
41088,
0,
'LMT',
    ],
    [
59958189312, #    utc_start 1900-12-31 12:35:12 (Mon)
62127694800, #      utc_end 1969-09-30 13:00:00 (Tue)
59958228912, #  local_start 1900-12-31 23:35:12 (Mon)
62127734400, #    local_end 1969-10-01 00:00:00 (Wed)
39600,
0,
'MHT',
    ],
    [
62127694800, #    utc_start 1969-09-30 13:00:00 (Tue)
DateTime::TimeZone::INFINITY, #      utc_end
62127738000, #  local_start 1969-10-01 01:00:00 (Wed)
DateTime::TimeZone::INFINITY, #    local_end
43200,
0,
'MHT',
    ],
];

sub olson_version {'2015d'}

sub has_dst_changes {0}

sub _max_year {2025}

sub _new_instance {
    return shift->_init( @_, spans => $spans );
}



1;

