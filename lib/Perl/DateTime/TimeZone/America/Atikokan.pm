# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/8FT049ktOU/northamerica.  Olson data version 2015d
#
# Do not edit this file directly.
#
package DateTime::TimeZone::America::Atikokan;
$DateTime::TimeZone::America::Atikokan::VERSION = '1.90';
use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::America::Atikokan::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
59768949988, #      utc_end 1895-01-01 06:06:28 (Tue)
DateTime::TimeZone::NEG_INFINITY, #  local_start
59768928000, #    local_end 1895-01-01 00:00:00 (Tue)
-21988,
0,
'LMT',
    ],
    [
59768949988, #    utc_start 1895-01-01 06:06:28 (Tue)
60503616000, #      utc_end 1918-04-14 08:00:00 (Sun)
59768928388, #  local_start 1895-01-01 00:06:28 (Tue)
60503594400, #    local_end 1918-04-14 02:00:00 (Sun)
-21600,
0,
'CST',
    ],
    [
60503616000, #    utc_start 1918-04-14 08:00:00 (Sun)
60520546800, #      utc_end 1918-10-27 07:00:00 (Sun)
60503598000, #  local_start 1918-04-14 03:00:00 (Sun)
60520528800, #    local_end 1918-10-27 02:00:00 (Sun)
-18000,
1,
'CDT',
    ],
    [
60520546800, #    utc_start 1918-10-27 07:00:00 (Sun)
61212434400, #      utc_end 1940-09-29 06:00:00 (Sun)
60520525200, #  local_start 1918-10-27 01:00:00 (Sun)
61212412800, #    local_end 1940-09-29 00:00:00 (Sun)
-21600,
0,
'CST',
    ],
    [
61212434400, #    utc_start 1940-09-29 06:00:00 (Sun)
61255468800, #      utc_end 1942-02-09 08:00:00 (Mon)
61212416400, #  local_start 1940-09-29 01:00:00 (Sun)
61255450800, #    local_end 1942-02-09 03:00:00 (Mon)
-18000,
1,
'CDT',
    ],
    [
61255468800, #    utc_start 1942-02-09 08:00:00 (Mon)
61366287600, #      utc_end 1945-08-14 23:00:00 (Tue)
61255450800, #  local_start 1942-02-09 03:00:00 (Mon)
61366269600, #    local_end 1945-08-14 18:00:00 (Tue)
-18000,
1,
'CWT',
    ],
    [
61366287600, #    utc_start 1945-08-14 23:00:00 (Tue)
61370290800, #      utc_end 1945-09-30 07:00:00 (Sun)
61366269600, #  local_start 1945-08-14 18:00:00 (Tue)
61370272800, #    local_end 1945-09-30 02:00:00 (Sun)
-18000,
1,
'CPT',
    ],
    [
61370290800, #    utc_start 1945-09-30 07:00:00 (Sun)
DateTime::TimeZone::INFINITY, #      utc_end
61370272800, #  local_start 1945-09-30 02:00:00 (Sun)
DateTime::TimeZone::INFINITY, #    local_end
-18000,
0,
'EST',
    ],
];

sub olson_version {'2015d'}

sub has_dst_changes {4}

sub _max_year {2025}

sub _new_instance {
    return shift->_init( @_, spans => $spans );
}



1;
