# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/8FT049ktOU/asia.  Olson data version 2015d
#
# Do not edit this file directly.
#
package DateTime::TimeZone::Asia::Baku;
$DateTime::TimeZone::Asia::Baku::VERSION = '1.90';
use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::Asia::Baku::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
60694519236, #      utc_end 1924-05-01 20:40:36 (Thu)
DateTime::TimeZone::NEG_INFINITY, #  local_start
60694531200, #    local_end 1924-05-02 00:00:00 (Fri)
11964,
0,
'LMT',
    ],
    [
60694519236, #    utc_start 1924-05-01 20:40:36 (Thu)
61730542800, #      utc_end 1957-02-28 21:00:00 (Thu)
60694530036, #  local_start 1924-05-01 23:40:36 (Thu)
61730553600, #    local_end 1957-03-01 00:00:00 (Fri)
10800,
0,
'BAKT',
    ],
    [
61730542800, #    utc_start 1957-02-28 21:00:00 (Thu)
62490600000, #      utc_end 1981-03-31 20:00:00 (Tue)
61730557200, #  local_start 1957-03-01 01:00:00 (Fri)
62490614400, #    local_end 1981-04-01 00:00:00 (Wed)
14400,
0,
'BAKT',
    ],
    [
62490600000, #    utc_start 1981-03-31 20:00:00 (Tue)
62506407600, #      utc_end 1981-09-30 19:00:00 (Wed)
62490618000, #  local_start 1981-04-01 01:00:00 (Wed)
62506425600, #    local_end 1981-10-01 00:00:00 (Thu)
18000,
1,
'BAKST',
    ],
    [
62506407600, #    utc_start 1981-09-30 19:00:00 (Wed)
62522136000, #      utc_end 1982-03-31 20:00:00 (Wed)
62506422000, #  local_start 1981-09-30 23:00:00 (Wed)
62522150400, #    local_end 1982-04-01 00:00:00 (Thu)
14400,
0,
'BAKT',
    ],
    [
62522136000, #    utc_start 1982-03-31 20:00:00 (Wed)
62537943600, #      utc_end 1982-09-30 19:00:00 (Thu)
62522154000, #  local_start 1982-04-01 01:00:00 (Thu)
62537961600, #    local_end 1982-10-01 00:00:00 (Fri)
18000,
1,
'BAKST',
    ],
    [
62537943600, #    utc_start 1982-09-30 19:00:00 (Thu)
62553672000, #      utc_end 1983-03-31 20:00:00 (Thu)
62537958000, #  local_start 1982-09-30 23:00:00 (Thu)
62553686400, #    local_end 1983-04-01 00:00:00 (Fri)
14400,
0,
'BAKT',
    ],
    [
62553672000, #    utc_start 1983-03-31 20:00:00 (Thu)
62569479600, #      utc_end 1983-09-30 19:00:00 (Fri)
62553690000, #  local_start 1983-04-01 01:00:00 (Fri)
62569497600, #    local_end 1983-10-01 00:00:00 (Sat)
18000,
1,
'BAKST',
    ],
    [
62569479600, #    utc_start 1983-09-30 19:00:00 (Fri)
62585294400, #      utc_end 1984-03-31 20:00:00 (Sat)
62569494000, #  local_start 1983-09-30 23:00:00 (Fri)
62585308800, #    local_end 1984-04-01 00:00:00 (Sun)
14400,
0,
'BAKT',
    ],
    [
62585294400, #    utc_start 1984-03-31 20:00:00 (Sat)
62601026400, #      utc_end 1984-09-29 22:00:00 (Sat)
62585312400, #  local_start 1984-04-01 01:00:00 (Sun)
62601044400, #    local_end 1984-09-30 03:00:00 (Sun)
18000,
1,
'BAKST',
    ],
    [
62601026400, #    utc_start 1984-09-29 22:00:00 (Sat)
62616751200, #      utc_end 1985-03-30 22:00:00 (Sat)
62601040800, #  local_start 1984-09-30 02:00:00 (Sun)
62616765600, #    local_end 1985-03-31 02:00:00 (Sun)
14400,
0,
'BAKT',
    ],
    [
62616751200, #    utc_start 1985-03-30 22:00:00 (Sat)
62632476000, #      utc_end 1985-09-28 22:00:00 (Sat)
62616769200, #  local_start 1985-03-31 03:00:00 (Sun)
62632494000, #    local_end 1985-09-29 03:00:00 (Sun)
18000,
1,
'BAKST',
    ],
    [
62632476000, #    utc_start 1985-09-28 22:00:00 (Sat)
62648200800, #      utc_end 1986-03-29 22:00:00 (Sat)
62632490400, #  local_start 1985-09-29 02:00:00 (Sun)
62648215200, #    local_end 1986-03-30 02:00:00 (Sun)
14400,
0,
'BAKT',
    ],
    [
62648200800, #    utc_start 1986-03-29 22:00:00 (Sat)
62663925600, #      utc_end 1986-09-27 22:00:00 (Sat)
62648218800, #  local_start 1986-03-30 03:00:00 (Sun)
62663943600, #    local_end 1986-09-28 03:00:00 (Sun)
18000,
1,
'BAKST',
    ],
    [
62663925600, #    utc_start 1986-09-27 22:00:00 (Sat)
62679650400, #      utc_end 1987-03-28 22:00:00 (Sat)
62663940000, #  local_start 1986-09-28 02:00:00 (Sun)
62679664800, #    local_end 1987-03-29 02:00:00 (Sun)
14400,
0,
'BAKT',
    ],
    [
62679650400, #    utc_start 1987-03-28 22:00:00 (Sat)
62695375200, #      utc_end 1987-09-26 22:00:00 (Sat)
62679668400, #  local_start 1987-03-29 03:00:00 (Sun)
62695393200, #    local_end 1987-09-27 03:00:00 (Sun)
18000,
1,
'BAKST',
    ],
    [
62695375200, #    utc_start 1987-09-26 22:00:00 (Sat)
62711100000, #      utc_end 1988-03-26 22:00:00 (Sat)
62695389600, #  local_start 1987-09-27 02:00:00 (Sun)
62711114400, #    local_end 1988-03-27 02:00:00 (Sun)
14400,
0,
'BAKT',
    ],
    [
62711100000, #    utc_start 1988-03-26 22:00:00 (Sat)
62726824800, #      utc_end 1988-09-24 22:00:00 (Sat)
62711118000, #  local_start 1988-03-27 03:00:00 (Sun)
62726842800, #    local_end 1988-09-25 03:00:00 (Sun)
18000,
1,
'BAKST',
    ],
    [
62726824800, #    utc_start 1988-09-24 22:00:00 (Sat)
62742549600, #      utc_end 1989-03-25 22:00:00 (Sat)
62726839200, #  local_start 1988-09-25 02:00:00 (Sun)
62742564000, #    local_end 1989-03-26 02:00:00 (Sun)
14400,
0,
'BAKT',
    ],
    [
62742549600, #    utc_start 1989-03-25 22:00:00 (Sat)
62758274400, #      utc_end 1989-09-23 22:00:00 (Sat)
62742567600, #  local_start 1989-03-26 03:00:00 (Sun)
62758292400, #    local_end 1989-09-24 03:00:00 (Sun)
18000,
1,
'BAKST',
    ],
    [
62758274400, #    utc_start 1989-09-23 22:00:00 (Sat)
62773999200, #      utc_end 1990-03-24 22:00:00 (Sat)
62758288800, #  local_start 1989-09-24 02:00:00 (Sun)
62774013600, #    local_end 1990-03-25 02:00:00 (Sun)
14400,
0,
'BAKT',
    ],
    [
62773999200, #    utc_start 1990-03-24 22:00:00 (Sat)
62790328800, #      utc_end 1990-09-29 22:00:00 (Sat)
62774017200, #  local_start 1990-03-25 03:00:00 (Sun)
62790346800, #    local_end 1990-09-30 03:00:00 (Sun)
18000,
1,
'BAKST',
    ],
    [
62790328800, #    utc_start 1990-09-29 22:00:00 (Sat)
62806053600, #      utc_end 1991-03-30 22:00:00 (Sat)
62790343200, #  local_start 1990-09-30 02:00:00 (Sun)
62806068000, #    local_end 1991-03-31 02:00:00 (Sun)
14400,
0,
'BAKT',
    ],
    [
62806053600, #    utc_start 1991-03-30 22:00:00 (Sat)
62819179200, #      utc_end 1991-08-29 20:00:00 (Thu)
62806068000, #  local_start 1991-03-31 02:00:00 (Sun)
62819193600, #    local_end 1991-08-30 00:00:00 (Fri)
14400,
1,
'BAKST',
    ],
    [
62819179200, #    utc_start 1991-08-29 20:00:00 (Thu)
62821782000, #      utc_end 1991-09-28 23:00:00 (Sat)
62819193600, #  local_start 1991-08-30 00:00:00 (Fri)
62821796400, #    local_end 1991-09-29 03:00:00 (Sun)
14400,
1,
'AZST',
    ],
    [
62821782000, #    utc_start 1991-09-28 23:00:00 (Sat)
62837496000, #      utc_end 1992-03-28 20:00:00 (Sat)
62821792800, #  local_start 1991-09-29 02:00:00 (Sun)
62837506800, #    local_end 1992-03-28 23:00:00 (Sat)
10800,
0,
'AZT',
    ],
    [
62837496000, #    utc_start 1992-03-28 20:00:00 (Sat)
62853217200, #      utc_end 1992-09-26 19:00:00 (Sat)
62837510400, #  local_start 1992-03-29 00:00:00 (Sun)
62853231600, #    local_end 1992-09-26 23:00:00 (Sat)
14400,
1,
'AZST',
    ],
    [
62853217200, #    utc_start 1992-09-26 19:00:00 (Sat)
62956123200, #      utc_end 1995-12-31 20:00:00 (Sun)
62853231600, #  local_start 1992-09-26 23:00:00 (Sat)
62956137600, #    local_end 1996-01-01 00:00:00 (Mon)
14400,
0,
'AZT',
    ],
    [
62956123200, #    utc_start 1995-12-31 20:00:00 (Sun)
62963917200, #      utc_end 1996-03-31 01:00:00 (Sun)
62956137600, #  local_start 1996-01-01 00:00:00 (Mon)
62963931600, #    local_end 1996-03-31 05:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
62963917200, #    utc_start 1996-03-31 01:00:00 (Sun)
62982061200, #      utc_end 1996-10-27 01:00:00 (Sun)
62963935200, #  local_start 1996-03-31 06:00:00 (Sun)
62982079200, #    local_end 1996-10-27 06:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
62982061200, #    utc_start 1996-10-27 01:00:00 (Sun)
62987745600, #      utc_end 1996-12-31 20:00:00 (Tue)
62982075600, #  local_start 1996-10-27 05:00:00 (Sun)
62987760000, #    local_end 1997-01-01 00:00:00 (Wed)
14400,
0,
'AZT',
    ],
    [
62987745600, #    utc_start 1996-12-31 20:00:00 (Tue)
62995363200, #      utc_end 1997-03-30 00:00:00 (Sun)
62987760000, #  local_start 1997-01-01 00:00:00 (Wed)
62995377600, #    local_end 1997-03-30 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
62995363200, #    utc_start 1997-03-30 00:00:00 (Sun)
63013507200, #      utc_end 1997-10-26 00:00:00 (Sun)
62995381200, #  local_start 1997-03-30 05:00:00 (Sun)
63013525200, #    local_end 1997-10-26 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63013507200, #    utc_start 1997-10-26 00:00:00 (Sun)
63026812800, #      utc_end 1998-03-29 00:00:00 (Sun)
63013521600, #  local_start 1997-10-26 04:00:00 (Sun)
63026827200, #    local_end 1998-03-29 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63026812800, #    utc_start 1998-03-29 00:00:00 (Sun)
63044956800, #      utc_end 1998-10-25 00:00:00 (Sun)
63026830800, #  local_start 1998-03-29 05:00:00 (Sun)
63044974800, #    local_end 1998-10-25 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63044956800, #    utc_start 1998-10-25 00:00:00 (Sun)
63058262400, #      utc_end 1999-03-28 00:00:00 (Sun)
63044971200, #  local_start 1998-10-25 04:00:00 (Sun)
63058276800, #    local_end 1999-03-28 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63058262400, #    utc_start 1999-03-28 00:00:00 (Sun)
63077011200, #      utc_end 1999-10-31 00:00:00 (Sun)
63058280400, #  local_start 1999-03-28 05:00:00 (Sun)
63077029200, #    local_end 1999-10-31 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63077011200, #    utc_start 1999-10-31 00:00:00 (Sun)
63089712000, #      utc_end 2000-03-26 00:00:00 (Sun)
63077025600, #  local_start 1999-10-31 04:00:00 (Sun)
63089726400, #    local_end 2000-03-26 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63089712000, #    utc_start 2000-03-26 00:00:00 (Sun)
63108460800, #      utc_end 2000-10-29 00:00:00 (Sun)
63089730000, #  local_start 2000-03-26 05:00:00 (Sun)
63108478800, #    local_end 2000-10-29 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63108460800, #    utc_start 2000-10-29 00:00:00 (Sun)
63121161600, #      utc_end 2001-03-25 00:00:00 (Sun)
63108475200, #  local_start 2000-10-29 04:00:00 (Sun)
63121176000, #    local_end 2001-03-25 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63121161600, #    utc_start 2001-03-25 00:00:00 (Sun)
63139910400, #      utc_end 2001-10-28 00:00:00 (Sun)
63121179600, #  local_start 2001-03-25 05:00:00 (Sun)
63139928400, #    local_end 2001-10-28 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63139910400, #    utc_start 2001-10-28 00:00:00 (Sun)
63153216000, #      utc_end 2002-03-31 00:00:00 (Sun)
63139924800, #  local_start 2001-10-28 04:00:00 (Sun)
63153230400, #    local_end 2002-03-31 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63153216000, #    utc_start 2002-03-31 00:00:00 (Sun)
63171360000, #      utc_end 2002-10-27 00:00:00 (Sun)
63153234000, #  local_start 2002-03-31 05:00:00 (Sun)
63171378000, #    local_end 2002-10-27 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63171360000, #    utc_start 2002-10-27 00:00:00 (Sun)
63184665600, #      utc_end 2003-03-30 00:00:00 (Sun)
63171374400, #  local_start 2002-10-27 04:00:00 (Sun)
63184680000, #    local_end 2003-03-30 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63184665600, #    utc_start 2003-03-30 00:00:00 (Sun)
63202809600, #      utc_end 2003-10-26 00:00:00 (Sun)
63184683600, #  local_start 2003-03-30 05:00:00 (Sun)
63202827600, #    local_end 2003-10-26 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63202809600, #    utc_start 2003-10-26 00:00:00 (Sun)
63216115200, #      utc_end 2004-03-28 00:00:00 (Sun)
63202824000, #  local_start 2003-10-26 04:00:00 (Sun)
63216129600, #    local_end 2004-03-28 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63216115200, #    utc_start 2004-03-28 00:00:00 (Sun)
63234864000, #      utc_end 2004-10-31 00:00:00 (Sun)
63216133200, #  local_start 2004-03-28 05:00:00 (Sun)
63234882000, #    local_end 2004-10-31 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63234864000, #    utc_start 2004-10-31 00:00:00 (Sun)
63247564800, #      utc_end 2005-03-27 00:00:00 (Sun)
63234878400, #  local_start 2004-10-31 04:00:00 (Sun)
63247579200, #    local_end 2005-03-27 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63247564800, #    utc_start 2005-03-27 00:00:00 (Sun)
63266313600, #      utc_end 2005-10-30 00:00:00 (Sun)
63247582800, #  local_start 2005-03-27 05:00:00 (Sun)
63266331600, #    local_end 2005-10-30 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63266313600, #    utc_start 2005-10-30 00:00:00 (Sun)
63279014400, #      utc_end 2006-03-26 00:00:00 (Sun)
63266328000, #  local_start 2005-10-30 04:00:00 (Sun)
63279028800, #    local_end 2006-03-26 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63279014400, #    utc_start 2006-03-26 00:00:00 (Sun)
63297763200, #      utc_end 2006-10-29 00:00:00 (Sun)
63279032400, #  local_start 2006-03-26 05:00:00 (Sun)
63297781200, #    local_end 2006-10-29 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63297763200, #    utc_start 2006-10-29 00:00:00 (Sun)
63310464000, #      utc_end 2007-03-25 00:00:00 (Sun)
63297777600, #  local_start 2006-10-29 04:00:00 (Sun)
63310478400, #    local_end 2007-03-25 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63310464000, #    utc_start 2007-03-25 00:00:00 (Sun)
63329212800, #      utc_end 2007-10-28 00:00:00 (Sun)
63310482000, #  local_start 2007-03-25 05:00:00 (Sun)
63329230800, #    local_end 2007-10-28 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63329212800, #    utc_start 2007-10-28 00:00:00 (Sun)
63342518400, #      utc_end 2008-03-30 00:00:00 (Sun)
63329227200, #  local_start 2007-10-28 04:00:00 (Sun)
63342532800, #    local_end 2008-03-30 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63342518400, #    utc_start 2008-03-30 00:00:00 (Sun)
63360662400, #      utc_end 2008-10-26 00:00:00 (Sun)
63342536400, #  local_start 2008-03-30 05:00:00 (Sun)
63360680400, #    local_end 2008-10-26 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63360662400, #    utc_start 2008-10-26 00:00:00 (Sun)
63373968000, #      utc_end 2009-03-29 00:00:00 (Sun)
63360676800, #  local_start 2008-10-26 04:00:00 (Sun)
63373982400, #    local_end 2009-03-29 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63373968000, #    utc_start 2009-03-29 00:00:00 (Sun)
63392112000, #      utc_end 2009-10-25 00:00:00 (Sun)
63373986000, #  local_start 2009-03-29 05:00:00 (Sun)
63392130000, #    local_end 2009-10-25 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63392112000, #    utc_start 2009-10-25 00:00:00 (Sun)
63405417600, #      utc_end 2010-03-28 00:00:00 (Sun)
63392126400, #  local_start 2009-10-25 04:00:00 (Sun)
63405432000, #    local_end 2010-03-28 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63405417600, #    utc_start 2010-03-28 00:00:00 (Sun)
63424166400, #      utc_end 2010-10-31 00:00:00 (Sun)
63405435600, #  local_start 2010-03-28 05:00:00 (Sun)
63424184400, #    local_end 2010-10-31 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63424166400, #    utc_start 2010-10-31 00:00:00 (Sun)
63436867200, #      utc_end 2011-03-27 00:00:00 (Sun)
63424180800, #  local_start 2010-10-31 04:00:00 (Sun)
63436881600, #    local_end 2011-03-27 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63436867200, #    utc_start 2011-03-27 00:00:00 (Sun)
63455616000, #      utc_end 2011-10-30 00:00:00 (Sun)
63436885200, #  local_start 2011-03-27 05:00:00 (Sun)
63455634000, #    local_end 2011-10-30 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63455616000, #    utc_start 2011-10-30 00:00:00 (Sun)
63468316800, #      utc_end 2012-03-25 00:00:00 (Sun)
63455630400, #  local_start 2011-10-30 04:00:00 (Sun)
63468331200, #    local_end 2012-03-25 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63468316800, #    utc_start 2012-03-25 00:00:00 (Sun)
63487065600, #      utc_end 2012-10-28 00:00:00 (Sun)
63468334800, #  local_start 2012-03-25 05:00:00 (Sun)
63487083600, #    local_end 2012-10-28 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63487065600, #    utc_start 2012-10-28 00:00:00 (Sun)
63500371200, #      utc_end 2013-03-31 00:00:00 (Sun)
63487080000, #  local_start 2012-10-28 04:00:00 (Sun)
63500385600, #    local_end 2013-03-31 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63500371200, #    utc_start 2013-03-31 00:00:00 (Sun)
63518515200, #      utc_end 2013-10-27 00:00:00 (Sun)
63500389200, #  local_start 2013-03-31 05:00:00 (Sun)
63518533200, #    local_end 2013-10-27 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63518515200, #    utc_start 2013-10-27 00:00:00 (Sun)
63531820800, #      utc_end 2014-03-30 00:00:00 (Sun)
63518529600, #  local_start 2013-10-27 04:00:00 (Sun)
63531835200, #    local_end 2014-03-30 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63531820800, #    utc_start 2014-03-30 00:00:00 (Sun)
63549964800, #      utc_end 2014-10-26 00:00:00 (Sun)
63531838800, #  local_start 2014-03-30 05:00:00 (Sun)
63549982800, #    local_end 2014-10-26 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63549964800, #    utc_start 2014-10-26 00:00:00 (Sun)
63563270400, #      utc_end 2015-03-29 00:00:00 (Sun)
63549979200, #  local_start 2014-10-26 04:00:00 (Sun)
63563284800, #    local_end 2015-03-29 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63563270400, #    utc_start 2015-03-29 00:00:00 (Sun)
63581414400, #      utc_end 2015-10-25 00:00:00 (Sun)
63563288400, #  local_start 2015-03-29 05:00:00 (Sun)
63581432400, #    local_end 2015-10-25 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63581414400, #    utc_start 2015-10-25 00:00:00 (Sun)
63594720000, #      utc_end 2016-03-27 00:00:00 (Sun)
63581428800, #  local_start 2015-10-25 04:00:00 (Sun)
63594734400, #    local_end 2016-03-27 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63594720000, #    utc_start 2016-03-27 00:00:00 (Sun)
63613468800, #      utc_end 2016-10-30 00:00:00 (Sun)
63594738000, #  local_start 2016-03-27 05:00:00 (Sun)
63613486800, #    local_end 2016-10-30 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63613468800, #    utc_start 2016-10-30 00:00:00 (Sun)
63626169600, #      utc_end 2017-03-26 00:00:00 (Sun)
63613483200, #  local_start 2016-10-30 04:00:00 (Sun)
63626184000, #    local_end 2017-03-26 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63626169600, #    utc_start 2017-03-26 00:00:00 (Sun)
63644918400, #      utc_end 2017-10-29 00:00:00 (Sun)
63626187600, #  local_start 2017-03-26 05:00:00 (Sun)
63644936400, #    local_end 2017-10-29 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63644918400, #    utc_start 2017-10-29 00:00:00 (Sun)
63657619200, #      utc_end 2018-03-25 00:00:00 (Sun)
63644932800, #  local_start 2017-10-29 04:00:00 (Sun)
63657633600, #    local_end 2018-03-25 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63657619200, #    utc_start 2018-03-25 00:00:00 (Sun)
63676368000, #      utc_end 2018-10-28 00:00:00 (Sun)
63657637200, #  local_start 2018-03-25 05:00:00 (Sun)
63676386000, #    local_end 2018-10-28 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63676368000, #    utc_start 2018-10-28 00:00:00 (Sun)
63689673600, #      utc_end 2019-03-31 00:00:00 (Sun)
63676382400, #  local_start 2018-10-28 04:00:00 (Sun)
63689688000, #    local_end 2019-03-31 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63689673600, #    utc_start 2019-03-31 00:00:00 (Sun)
63707817600, #      utc_end 2019-10-27 00:00:00 (Sun)
63689691600, #  local_start 2019-03-31 05:00:00 (Sun)
63707835600, #    local_end 2019-10-27 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63707817600, #    utc_start 2019-10-27 00:00:00 (Sun)
63721123200, #      utc_end 2020-03-29 00:00:00 (Sun)
63707832000, #  local_start 2019-10-27 04:00:00 (Sun)
63721137600, #    local_end 2020-03-29 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63721123200, #    utc_start 2020-03-29 00:00:00 (Sun)
63739267200, #      utc_end 2020-10-25 00:00:00 (Sun)
63721141200, #  local_start 2020-03-29 05:00:00 (Sun)
63739285200, #    local_end 2020-10-25 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63739267200, #    utc_start 2020-10-25 00:00:00 (Sun)
63752572800, #      utc_end 2021-03-28 00:00:00 (Sun)
63739281600, #  local_start 2020-10-25 04:00:00 (Sun)
63752587200, #    local_end 2021-03-28 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63752572800, #    utc_start 2021-03-28 00:00:00 (Sun)
63771321600, #      utc_end 2021-10-31 00:00:00 (Sun)
63752590800, #  local_start 2021-03-28 05:00:00 (Sun)
63771339600, #    local_end 2021-10-31 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63771321600, #    utc_start 2021-10-31 00:00:00 (Sun)
63784022400, #      utc_end 2022-03-27 00:00:00 (Sun)
63771336000, #  local_start 2021-10-31 04:00:00 (Sun)
63784036800, #    local_end 2022-03-27 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63784022400, #    utc_start 2022-03-27 00:00:00 (Sun)
63802771200, #      utc_end 2022-10-30 00:00:00 (Sun)
63784040400, #  local_start 2022-03-27 05:00:00 (Sun)
63802789200, #    local_end 2022-10-30 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63802771200, #    utc_start 2022-10-30 00:00:00 (Sun)
63815472000, #      utc_end 2023-03-26 00:00:00 (Sun)
63802785600, #  local_start 2022-10-30 04:00:00 (Sun)
63815486400, #    local_end 2023-03-26 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63815472000, #    utc_start 2023-03-26 00:00:00 (Sun)
63834220800, #      utc_end 2023-10-29 00:00:00 (Sun)
63815490000, #  local_start 2023-03-26 05:00:00 (Sun)
63834238800, #    local_end 2023-10-29 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63834220800, #    utc_start 2023-10-29 00:00:00 (Sun)
63847526400, #      utc_end 2024-03-31 00:00:00 (Sun)
63834235200, #  local_start 2023-10-29 04:00:00 (Sun)
63847540800, #    local_end 2024-03-31 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63847526400, #    utc_start 2024-03-31 00:00:00 (Sun)
63865670400, #      utc_end 2024-10-27 00:00:00 (Sun)
63847544400, #  local_start 2024-03-31 05:00:00 (Sun)
63865688400, #    local_end 2024-10-27 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63865670400, #    utc_start 2024-10-27 00:00:00 (Sun)
63878976000, #      utc_end 2025-03-30 00:00:00 (Sun)
63865684800, #  local_start 2024-10-27 04:00:00 (Sun)
63878990400, #    local_end 2025-03-30 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63878976000, #    utc_start 2025-03-30 00:00:00 (Sun)
63897120000, #      utc_end 2025-10-26 00:00:00 (Sun)
63878994000, #  local_start 2025-03-30 05:00:00 (Sun)
63897138000, #    local_end 2025-10-26 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
    [
63897120000, #    utc_start 2025-10-26 00:00:00 (Sun)
63910425600, #      utc_end 2026-03-29 00:00:00 (Sun)
63897134400, #  local_start 2025-10-26 04:00:00 (Sun)
63910440000, #    local_end 2026-03-29 04:00:00 (Sun)
14400,
0,
'AZT',
    ],
    [
63910425600, #    utc_start 2026-03-29 00:00:00 (Sun)
63928569600, #      utc_end 2026-10-25 00:00:00 (Sun)
63910443600, #  local_start 2026-03-29 05:00:00 (Sun)
63928587600, #    local_end 2026-10-25 05:00:00 (Sun)
18000,
1,
'AZST',
    ],
];

sub olson_version {'2015d'}

sub has_dst_changes {44}

sub _max_year {2025}

sub _new_instance {
    return shift->_init( @_, spans => $spans );
}

sub _last_offset { 14400 }

my $last_observance = bless( {
  'format' => 'AZ%sT',
  'gmtoff' => '4:00',
  'local_start_datetime' => bless( {
    'formatter' => undef,
    'local_rd_days' => 729025,
    'local_rd_secs' => 0,
    'offset_modifier' => 0,
    'rd_nanosecs' => 0,
    'tz' => bless( {
      'name' => 'floating',
      'offset' => 0
    }, 'DateTime::TimeZone::Floating' ),
    'utc_rd_days' => 729025,
    'utc_rd_secs' => 0,
    'utc_year' => 1998
  }, 'DateTime' ),
  'offset_from_std' => 0,
  'offset_from_utc' => 14400,
  'until' => [],
  'utc_start_datetime' => bless( {
    'formatter' => undef,
    'local_rd_days' => 729024,
    'local_rd_secs' => 72000,
    'offset_modifier' => 0,
    'rd_nanosecs' => 0,
    'tz' => bless( {
      'name' => 'floating',
      'offset' => 0
    }, 'DateTime::TimeZone::Floating' ),
    'utc_rd_days' => 729024,
    'utc_rd_secs' => 72000,
    'utc_year' => 1997
  }, 'DateTime' )
}, 'DateTime::TimeZone::OlsonDB::Observance' )
;
sub _last_observance { $last_observance }

my $rules = [
  bless( {
    'at' => '4:00',
    'from' => '1997',
    'in' => 'Mar',
    'letter' => 'S',
    'name' => 'Azer',
    'offset_from_std' => 3600,
    'on' => 'lastSun',
    'save' => '1:00',
    'to' => 'max',
    'type' => undef
  }, 'DateTime::TimeZone::OlsonDB::Rule' ),
  bless( {
    'at' => '5:00',
    'from' => '1997',
    'in' => 'Oct',
    'letter' => '',
    'name' => 'Azer',
    'offset_from_std' => 0,
    'on' => 'lastSun',
    'save' => '0',
    'to' => 'max',
    'type' => undef
  }, 'DateTime::TimeZone::OlsonDB::Rule' )
]
;
sub _rules { $rules }


1;

