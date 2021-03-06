###########################################################################
#
# This file is auto-generated by the Perl DateTime Suite locale
# generator (0.05).  This code generator comes with the
# DateTime::Locale distribution in the tools/ directory, and is called
# generate-from-cldr.
#
# This file as generated from the CLDR XML locale data.  See the
# LICENSE.cldr file included in this distribution for license details.
#
# This file was generated from the source file en_Shaw.xml
# The source file version number was 1.9, generated on
# 2009/05/05 23:06:35.
#
# Do not edit this file directly.
#
###########################################################################

package DateTime::Locale::en_Shaw;

use strict;
use warnings;

our $VERSION = '0.46';

use utf8;

use base 'DateTime::Locale::en';

sub cldr_version { return "1\.7\.1" }

{
    my $am_pm_abbreviated = [ "饜懆饜懃", "饜憪饜懃" ];
    sub am_pm_abbreviated { return $am_pm_abbreviated }
}
{
    my $day_format_abbreviated = [ "路饜懃饜懎", "路饜憫饜懙", "路饜憿饜懅", "路饜憯饜懟", "路饜憮饜懏", "路饜憰饜懆", "路饜憰饜懎" ];
    sub day_format_abbreviated { return $day_format_abbreviated }
}

sub day_format_narrow { $_[0]->day_stand_alone_narrow() }

{
    my $day_format_wide = [ "路饜懃饜懎饜憴饜憶饜懕", "路饜憫饜懙饜憻饜憶饜懕", "路饜憿饜懅饜憴饜憻饜憶饜懕", "路饜憯饜懟饜憻饜憶饜懕", "路饜憮饜懏饜懖饜憶饜懕", "路饜憰饜懆饜憶饜懟饜憶饜懕", "路饜憰饜懎饜憴饜憶饜懕" ];
    sub day_format_wide { return $day_format_wide }
}

sub day_stand_alone_abbreviated { $_[0]->day_format_abbreviated() }

{
    my $day_stand_alone_narrow = [ "饜懃", "饜憫", "饜憿", "饜憯", "饜憮", "饜憰", "饜憰" ];
    sub day_stand_alone_narrow { return $day_stand_alone_narrow }
}

sub day_stand_alone_wide { $_[0]->day_format_wide() }

{
    my $era_abbreviated = [ "饜憵路饜憭", "饜懆饜憶" ];
    sub era_abbreviated { return $era_abbreviated }
}
{
    my $era_narrow = [ "饜憵", "饜懆" ];
    sub era_narrow { return $era_narrow }
}
{
    my $era_wide = [ "饜憵饜懓饜憮饜應饜懏\ 路饜憭饜懏饜懖饜憰饜憫", "饜懆饜憴饜懘\ 饜憶饜應饜懃饜懄饜憴饜懓" ];
    sub era_wide { return $era_wide }
}
{
    my $first_day_of_week = "1";
    sub first_day_of_week { return $first_day_of_week }
}

{
    my $month_format_abbreviated = [ "路饜憽饜懆", "路饜憮饜懅", "路饜懃饜懜", "路饜懕饜憪", "路饜懃饜懕", "路饜憽饜懙", "路饜憽饜懌", "路饜應饜憸", "路饜憰饜懅", "路饜懛饜憭", "路饜懐饜懘", "路饜憶饜懎" ];
    sub month_format_abbreviated { return $month_format_abbreviated }
}

sub month_format_narrow { $_[0]->month_stand_alone_narrow() }

{
    my $month_format_wide = [ "路饜憽饜懆饜憴饜憳饜懎饜憿饜懞饜懓", "路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓", "路饜懃饜懜饜憲", "路饜懕饜憪饜懏饜懎饜懁", "路饜懃饜懕", "路饜憽饜懙饜懐", "路饜憽饜懌饜懁饜懖", "路饜應饜憸饜懎饜憰饜憫", "路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜", "路饜懛饜憭饜憫饜懘饜憵饜懜", "路饜懐饜懘饜憹饜懅饜懃饜憵饜懜", "路饜憶饜懎饜憰饜懅饜懃饜憵饜懜" ];
    sub month_format_wide { return $month_format_wide }
}

sub month_stand_alone_abbreviated { $_[0]->month_format_abbreviated() }

{
    my $month_stand_alone_narrow = [ "饜憽", "饜憮", "饜懃", "饜懕", "饜懃", "饜憽", "饜憽", "饜應", "饜憰", "饜懛", "饜懐", "饜憶" ];
    sub month_stand_alone_narrow { return $month_stand_alone_narrow }
}

sub month_stand_alone_wide { $_[0]->month_format_wide() }

{
    my $quarter_format_abbreviated = [ "饜憭1", "饜憭2", "饜憭3", "饜憭4" ];
    sub quarter_format_abbreviated { return $quarter_format_abbreviated }
}
{
    my $quarter_format_wide = [ "1饜憰饜憫\ 饜憭饜憿饜懜饜憶饜懜", "2饜懐饜憶\ 饜憭饜憿饜懜饜憶饜懜", "3饜懟饜憶\ 饜憭饜憿饜懜饜憶饜懜", "4饜懝饜憯\ 饜憭饜憿饜懜饜憶饜懜" ];
    sub quarter_format_wide { return $quarter_format_wide }
}

sub quarter_stand_alone_abbreviated { $_[0]->quarter_format_abbreviated() }


sub quarter_stand_alone_wide { $_[0]->quarter_format_wide() }

1;

__END__


=pod

=encoding utf8

=head1 NAME

DateTime::Locale::en_Shaw

=head1 SYNOPSIS

  use DateTime;

  my $dt = DateTime->now( locale => 'en_Shaw' );
  print $dt->month_name();

=head1 DESCRIPTION

This is the DateTime locale package for English Shavian.

=head1 DATA

This locale inherits from the L<DateTime::Locale::en> locale.

It contains the following data.

=head2 Days

=head3 Wide (format)

  路饜懃饜懎饜憴饜憶饜懕
  路饜憫饜懙饜憻饜憶饜懕
  路饜憿饜懅饜憴饜憻饜憶饜懕
  路饜憯饜懟饜憻饜憶饜懕
  路饜憮饜懏饜懖饜憶饜懕
  路饜憰饜懆饜憶饜懟饜憶饜懕
  路饜憰饜懎饜憴饜憶饜懕

=head3 Abbreviated (format)

  路饜懃饜懎
  路饜憫饜懙
  路饜憿饜懅
  路饜憯饜懟
  路饜憮饜懏
  路饜憰饜懆
  路饜憰饜懎

=head3 Narrow (format)

  饜懃
  饜憫
  饜憿
  饜憯
  饜憮
  饜憰
  饜憰

=head3 Wide (stand-alone)

  路饜懃饜懎饜憴饜憶饜懕
  路饜憫饜懙饜憻饜憶饜懕
  路饜憿饜懅饜憴饜憻饜憶饜懕
  路饜憯饜懟饜憻饜憶饜懕
  路饜憮饜懏饜懖饜憶饜懕
  路饜憰饜懆饜憶饜懟饜憶饜懕
  路饜憰饜懎饜憴饜憶饜懕

=head3 Abbreviated (stand-alone)

  路饜懃饜懎
  路饜憫饜懙
  路饜憿饜懅
  路饜憯饜懟
  路饜憮饜懏
  路饜憰饜懆
  路饜憰饜懎

=head3 Narrow (stand-alone)

  饜懃
  饜憫
  饜憿
  饜憯
  饜憮
  饜憰
  饜憰

=head2 Months

=head3 Wide (format)

  路饜憽饜懆饜憴饜憳饜懎饜憿饜懞饜懓
  路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓
  路饜懃饜懜饜憲
  路饜懕饜憪饜懏饜懎饜懁
  路饜懃饜懕
  路饜憽饜懙饜懐
  路饜憽饜懌饜懁饜懖
  路饜應饜憸饜懎饜憰饜憫
  路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜
  路饜懛饜憭饜憫饜懘饜憵饜懜
  路饜懐饜懘饜憹饜懅饜懃饜憵饜懜
  路饜憶饜懎饜憰饜懅饜懃饜憵饜懜

=head3 Abbreviated (format)

  路饜憽饜懆
  路饜憮饜懅
  路饜懃饜懜
  路饜懕饜憪
  路饜懃饜懕
  路饜憽饜懙
  路饜憽饜懌
  路饜應饜憸
  路饜憰饜懅
  路饜懛饜憭
  路饜懐饜懘
  路饜憶饜懎

=head3 Narrow (format)

  饜憽
  饜憮
  饜懃
  饜懕
  饜懃
  饜憽
  饜憽
  饜應
  饜憰
  饜懛
  饜懐
  饜憶

=head3 Wide (stand-alone)

  路饜憽饜懆饜憴饜憳饜懎饜憿饜懞饜懓
  路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓
  路饜懃饜懜饜憲
  路饜懕饜憪饜懏饜懎饜懁
  路饜懃饜懕
  路饜憽饜懙饜懐
  路饜憽饜懌饜懁饜懖
  路饜應饜憸饜懎饜憰饜憫
  路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜
  路饜懛饜憭饜憫饜懘饜憵饜懜
  路饜懐饜懘饜憹饜懅饜懃饜憵饜懜
  路饜憶饜懎饜憰饜懅饜懃饜憵饜懜

=head3 Abbreviated (stand-alone)

  路饜憽饜懆
  路饜憮饜懅
  路饜懃饜懜
  路饜懕饜憪
  路饜懃饜懕
  路饜憽饜懙
  路饜憽饜懌
  路饜應饜憸
  路饜憰饜懅
  路饜懛饜憭
  路饜懐饜懘
  路饜憶饜懎

=head3 Narrow (stand-alone)

  饜憽
  饜憮
  饜懃
  饜懕
  饜懃
  饜憽
  饜憽
  饜應
  饜憰
  饜懛
  饜懐
  饜憶

=head2 Quarters

=head3 Wide (format)

  1饜憰饜憫 饜憭饜憿饜懜饜憶饜懜
  2饜懐饜憶 饜憭饜憿饜懜饜憶饜懜
  3饜懟饜憶 饜憭饜憿饜懜饜憶饜懜
  4饜懝饜憯 饜憭饜憿饜懜饜憶饜懜

=head3 Abbreviated (format)

  饜憭1
  饜憭2
  饜憭3
  饜憭4

=head3 Narrow (format)

  1
  2
  3
  4

=head3 Wide (stand-alone)

  1饜憰饜憫 饜憭饜憿饜懜饜憶饜懜
  2饜懐饜憶 饜憭饜憿饜懜饜憶饜懜
  3饜懟饜憶 饜憭饜憿饜懜饜憶饜懜
  4饜懝饜憯 饜憭饜憿饜懜饜憶饜懜

=head3 Abbreviated (stand-alone)

  饜憭1
  饜憭2
  饜憭3
  饜憭4

=head3 Narrow (stand-alone)

  1
  2
  3
  4

=head2 Eras

=head3 Wide

  饜憵饜懓饜憮饜應饜懏 路饜憭饜懏饜懖饜憰饜憫
  饜懆饜憴饜懘 饜憶饜應饜懃饜懄饜憴饜懓

=head3 Abbreviated

  饜憵路饜憭
  饜懆饜憶

=head3 Narrow

  饜憵
  饜懆

=head2 Date Formats

=head3 Full

   2008-02-05T18:30:30 = 路饜憫饜懙饜憻饜憶饜懕, 路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓 5, 2008
   1995-12-22T09:05:02 = 路饜憮饜懏饜懖饜憶饜懕, 路饜憶饜懎饜憰饜懅饜懃饜憵饜懜 22, 1995
  -0010-09-15T04:44:23 = 路饜憰饜懆饜憶饜懟饜憶饜懕, 路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜 15, -10

=head3 Long

   2008-02-05T18:30:30 = 路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓 5, 2008
   1995-12-22T09:05:02 = 路饜憶饜懎饜憰饜懅饜懃饜憵饜懜 22, 1995
  -0010-09-15T04:44:23 = 路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜 15, -10

=head3 Medium

   2008-02-05T18:30:30 = 路饜憮饜懅 5, 2008
   1995-12-22T09:05:02 = 路饜憶饜懎 22, 1995
  -0010-09-15T04:44:23 = 路饜憰饜懅 15, -10

=head3 Short

   2008-02-05T18:30:30 = 2/5/08
   1995-12-22T09:05:02 = 12/22/95
  -0010-09-15T04:44:23 = 9/15/-10

=head3 Default

   2008-02-05T18:30:30 = 路饜憮饜懅 5, 2008
   1995-12-22T09:05:02 = 路饜憶饜懎 22, 1995
  -0010-09-15T04:44:23 = 路饜憰饜懅 15, -10

=head2 Time Formats

=head3 Full

   2008-02-05T18:30:30 = 6:30:30 饜憪饜懃 UTC
   1995-12-22T09:05:02 = 9:05:02 饜懆饜懃 UTC
  -0010-09-15T04:44:23 = 4:44:23 饜懆饜懃 UTC

=head3 Long

   2008-02-05T18:30:30 = 6:30:30 饜憪饜懃 UTC
   1995-12-22T09:05:02 = 9:05:02 饜懆饜懃 UTC
  -0010-09-15T04:44:23 = 4:44:23 饜懆饜懃 UTC

=head3 Medium

   2008-02-05T18:30:30 = 6:30:30 饜憪饜懃
   1995-12-22T09:05:02 = 9:05:02 饜懆饜懃
  -0010-09-15T04:44:23 = 4:44:23 饜懆饜懃

=head3 Short

   2008-02-05T18:30:30 = 6:30 饜憪饜懃
   1995-12-22T09:05:02 = 9:05 饜懆饜懃
  -0010-09-15T04:44:23 = 4:44 饜懆饜懃

=head3 Default

   2008-02-05T18:30:30 = 6:30:30 饜憪饜懃
   1995-12-22T09:05:02 = 9:05:02 饜懆饜懃
  -0010-09-15T04:44:23 = 4:44:23 饜懆饜懃

=head2 Datetime Formats

=head3 Full

   2008-02-05T18:30:30 = 路饜憫饜懙饜憻饜憶饜懕, 路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓 5, 2008 6:30:30 饜憪饜懃 UTC
   1995-12-22T09:05:02 = 路饜憮饜懏饜懖饜憶饜懕, 路饜憶饜懎饜憰饜懅饜懃饜憵饜懜 22, 1995 9:05:02 饜懆饜懃 UTC
  -0010-09-15T04:44:23 = 路饜憰饜懆饜憶饜懟饜憶饜懕, 路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜 15, -10 4:44:23 饜懆饜懃 UTC

=head3 Long

   2008-02-05T18:30:30 = 路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓 5, 2008 6:30:30 饜憪饜懃 UTC
   1995-12-22T09:05:02 = 路饜憶饜懎饜憰饜懅饜懃饜憵饜懜 22, 1995 9:05:02 饜懆饜懃 UTC
  -0010-09-15T04:44:23 = 路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜 15, -10 4:44:23 饜懆饜懃 UTC

=head3 Medium

   2008-02-05T18:30:30 = 路饜憮饜懅 5, 2008 6:30:30 饜憪饜懃
   1995-12-22T09:05:02 = 路饜憶饜懎 22, 1995 9:05:02 饜懆饜懃
  -0010-09-15T04:44:23 = 路饜憰饜懅 15, -10 4:44:23 饜懆饜懃

=head3 Short

   2008-02-05T18:30:30 = 2/5/08 6:30 饜憪饜懃
   1995-12-22T09:05:02 = 12/22/95 9:05 饜懆饜懃
  -0010-09-15T04:44:23 = 9/15/-10 4:44 饜懆饜懃

=head3 Default

   2008-02-05T18:30:30 = 路饜憮饜懅 5, 2008 6:30:30 饜憪饜懃
   1995-12-22T09:05:02 = 路饜憶饜懎 22, 1995 9:05:02 饜懆饜懃
  -0010-09-15T04:44:23 = 路饜憰饜懅 15, -10 4:44:23 饜懆饜懃

=head2 Available Formats

=head3 d (d)

   2008-02-05T18:30:30 = 5
   1995-12-22T09:05:02 = 22
  -0010-09-15T04:44:23 = 15

=head3 EEEd (d EEE)

   2008-02-05T18:30:30 = 5 路饜憫饜懙
   1995-12-22T09:05:02 = 22 路饜憮饜懏
  -0010-09-15T04:44:23 = 15 路饜憰饜懆

=head3 Hm (H:mm)

   2008-02-05T18:30:30 = 18:30
   1995-12-22T09:05:02 = 9:05
  -0010-09-15T04:44:23 = 4:44

=head3 hm (h:mm a)

   2008-02-05T18:30:30 = 6:30 饜憪饜懃
   1995-12-22T09:05:02 = 9:05 饜懆饜懃
  -0010-09-15T04:44:23 = 4:44 饜懆饜懃

=head3 Hms (H:mm:ss)

   2008-02-05T18:30:30 = 18:30:30
   1995-12-22T09:05:02 = 9:05:02
  -0010-09-15T04:44:23 = 4:44:23

=head3 hms (h:mm:ss a)

   2008-02-05T18:30:30 = 6:30:30 饜憪饜懃
   1995-12-22T09:05:02 = 9:05:02 饜懆饜懃
  -0010-09-15T04:44:23 = 4:44:23 饜懆饜懃

=head3 M (L)

   2008-02-05T18:30:30 = 2
   1995-12-22T09:05:02 = 12
  -0010-09-15T04:44:23 = 9

=head3 Md (M/d)

   2008-02-05T18:30:30 = 2/5
   1995-12-22T09:05:02 = 12/22
  -0010-09-15T04:44:23 = 9/15

=head3 MEd (E, M/d)

   2008-02-05T18:30:30 = 路饜憫饜懙, 2/5
   1995-12-22T09:05:02 = 路饜憮饜懏, 12/22
  -0010-09-15T04:44:23 = 路饜憰饜懆, 9/15

=head3 MMM (LLL)

   2008-02-05T18:30:30 = 路饜憮饜懅
   1995-12-22T09:05:02 = 路饜憶饜懎
  -0010-09-15T04:44:23 = 路饜憰饜懅

=head3 MMMd (MMM d)

   2008-02-05T18:30:30 = 路饜憮饜懅 5
   1995-12-22T09:05:02 = 路饜憶饜懎 22
  -0010-09-15T04:44:23 = 路饜憰饜懅 15

=head3 MMMEd (E, MMM d)

   2008-02-05T18:30:30 = 路饜憫饜懙, 路饜憮饜懅 5
   1995-12-22T09:05:02 = 路饜憮饜懏, 路饜憶饜懎 22
  -0010-09-15T04:44:23 = 路饜憰饜懆, 路饜憰饜懅 15

=head3 MMMMd (MMMM d)

   2008-02-05T18:30:30 = 路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓 5
   1995-12-22T09:05:02 = 路饜憶饜懎饜憰饜懅饜懃饜憵饜懜 22
  -0010-09-15T04:44:23 = 路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜 15

=head3 MMMMEd (E, MMMM d)

   2008-02-05T18:30:30 = 路饜憫饜懙, 路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓 5
   1995-12-22T09:05:02 = 路饜憮饜懏, 路饜憶饜懎饜憰饜懅饜懃饜憵饜懜 22
  -0010-09-15T04:44:23 = 路饜憰饜懆, 路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜 15

=head3 ms (mm:ss)

   2008-02-05T18:30:30 = 30:30
   1995-12-22T09:05:02 = 05:02
  -0010-09-15T04:44:23 = 44:23

=head3 y (y)

   2008-02-05T18:30:30 = 2008
   1995-12-22T09:05:02 = 1995
  -0010-09-15T04:44:23 = -10

=head3 yM (M/yyyy)

   2008-02-05T18:30:30 = 2/2008
   1995-12-22T09:05:02 = 12/1995
  -0010-09-15T04:44:23 = 9/-010

=head3 yMEd (EEE, M/d/yyyy)

   2008-02-05T18:30:30 = 路饜憫饜懙, 2/5/2008
   1995-12-22T09:05:02 = 路饜憮饜懏, 12/22/1995
  -0010-09-15T04:44:23 = 路饜憰饜懆, 9/15/-010

=head3 yMMM (MMM y)

   2008-02-05T18:30:30 = 路饜憮饜懅 2008
   1995-12-22T09:05:02 = 路饜憶饜懎 1995
  -0010-09-15T04:44:23 = 路饜憰饜懅 -10

=head3 yMMMEd (EEE, MMM d, y)

   2008-02-05T18:30:30 = 路饜憫饜懙, 路饜憮饜懅 5, 2008
   1995-12-22T09:05:02 = 路饜憮饜懏, 路饜憶饜懎 22, 1995
  -0010-09-15T04:44:23 = 路饜憰饜懆, 路饜憰饜懅 15, -10

=head3 yMMMM (MMMM y)

   2008-02-05T18:30:30 = 路饜憮饜懅饜憵饜憳饜懙饜憿饜懞饜懓 2008
   1995-12-22T09:05:02 = 路饜憶饜懎饜憰饜懅饜懃饜憵饜懜 1995
  -0010-09-15T04:44:23 = 路饜憰饜懅饜憪饜憫饜懅饜懃饜憵饜懜 -10

=head3 yQ (Q yyyy)

   2008-02-05T18:30:30 = 1 2008
   1995-12-22T09:05:02 = 4 1995
  -0010-09-15T04:44:23 = 3 -010

=head3 yQQQ (QQQ y)

   2008-02-05T18:30:30 = 饜憭1 2008
   1995-12-22T09:05:02 = 饜憭4 1995
  -0010-09-15T04:44:23 = 饜憭3 -10

=head2 Miscellaneous

=head3 Prefers 24 hour time?

No

=head3 Local first day of the week

路饜懃饜懎饜憴饜憶饜懕


=head1 SUPPORT

See L<DateTime::Locale>.

=cut
