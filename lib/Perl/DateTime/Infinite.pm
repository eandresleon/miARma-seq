package DateTime::Infinite;

use strict;
use warnings;

our $VERSION = '1.19';

use DateTime;
use DateTime::TimeZone;

use base qw(DateTime);

foreach my $m (qw( set set_time_zone truncate )) {
    no strict 'refs';
    *{"DateTime::Infinite::$m"} = sub { return $_[0] };
}

sub is_finite   {0}
sub is_infinite {1}

sub _rd2ymd {
    return $_[2] ? ( $_[1] ) x 7 : ( $_[1] ) x 3;
}

sub _seconds_as_components {
    return ( $_[1] ) x 3;
}

sub _stringify {
    $_[0]->{utc_rd_days} == DateTime::INFINITY
        ? DateTime::INFINITY . ''
        : DateTime::NEG_INFINITY . '';
}

sub STORABLE_freeze {return}
sub STORABLE_thaw   {return}

package DateTime::Infinite::Future;

use strict;
use warnings;

use base qw(DateTime::Infinite);

{
    my $Pos = bless {
        utc_rd_days   => DateTime::INFINITY,
        utc_rd_secs   => DateTime::INFINITY,
        local_rd_days => DateTime::INFINITY,
        local_rd_secs => DateTime::INFINITY,
        rd_nanosecs   => DateTime::INFINITY,
        tz            => DateTime::TimeZone->new( name => 'floating' ),
        locale        => FakeLocale->instance(),
        },
        __PACKAGE__;

    $Pos->_calc_utc_rd;
    $Pos->_calc_local_rd;

    sub new {$Pos}
}

package DateTime::Infinite::Past;

use strict;
use warnings;

use base qw(DateTime::Infinite);

{
    my $Neg = bless {
        utc_rd_days   => DateTime::NEG_INFINITY,
        utc_rd_secs   => DateTime::NEG_INFINITY,
        local_rd_days => DateTime::NEG_INFINITY,
        local_rd_secs => DateTime::NEG_INFINITY,
        rd_nanosecs   => DateTime::NEG_INFINITY,
        tz            => DateTime::TimeZone->new( name => 'floating' ),
        locale        => FakeLocale->instance(),
        },
        __PACKAGE__;

    $Neg->_calc_utc_rd;
    $Neg->_calc_local_rd;

    sub new {$Neg}
}

package    # hide from PAUSE
    FakeLocale;

use strict;
use warnings;

use DateTime::Locale;

my $Instance;

sub instance {
    return $Instance ||= bless { locale => DateTime::Locale->load('en_US') },
        __PACKAGE__;
}

sub id {
    return 'infinite';
}

sub language_id {
    return 'infinite';
}

sub name {
    'Fake locale for Infinite DateTime objects';
}

sub language {
    'Fake locale for Infinite DateTime objects';
}

my @methods = qw(
    script_id
    territory_id
    variant_id
    script
    territory
    variant
    native_name
    native_language
    native_script
    native_territory
    native_variant
);

for my $meth (@methods) {
    no strict 'refs';
    *{$meth} = sub {undef};
}

# Totally arbitrary
sub first_day_of_week {
    return 1;
}

sub prefers_24_hour_time {
    return 0;
}

our $AUTOLOAD;

sub AUTOLOAD {
    my $self = shift;

    my ($meth) = $AUTOLOAD =~ /::(\w+)$/;

    if ( $meth =~ /format/ && $meth !~ /^(?:day|month|quarter)/ ) {
        return $self->{locale}->$meth(@_);
    }

    return [];
}

1;

# ABSTRACT: Infinite past and future DateTime objects

__END__

=pod

=head1 NAME

DateTime::Infinite - Infinite past and future DateTime objects

=head1 VERSION

version 1.19

=head1 SYNOPSIS

  my $future = DateTime::Infinite::Future->new();
  my $past   = DateTime::Infinite::Past->new();

=head1 DESCRIPTION

This module provides two L<DateTime.pm|DateTime> subclasses,
C<DateTime::Infinite::Future> and C<DateTime::Infinite::Past>.

The objects are in the "floating" timezone, and this cannot be
changed.

=head1 BUGS

There seem to be lots of problems when dealing with infinite numbers
on Win32. This may be a problem with this code, Perl, or Win32's IEEE
math implementation. Either way, the module may not be well-behaved
on Win32 operating systems.

=head1 METHODS

The only constructor for these two classes is the C<new()> method, as
shown in the L<SYNOPSIS|/SYNOPSIS>. This method takes no parameters.

All "get" methods in this module simply return infinity, positive or
negative. If the method is expected to return a string, it return the
string representation of positive or negative infinity used by your
system. For example, on my system calling C<year()> returns a number
which when printed appears either "inf" or "-inf".

The object is not mutable, so the C<set()>, C<set_time_zone()>, and
C<truncate()> methods are all do-nothing methods that simply return
the object they are called with.

Obviously, the C<is_finite()> method returns false and the
C<is_infinite()> method returns true.

=head1 SEE ALSO

datetime@perl.org mailing list

http://datetime.perl.org/

=head1 AUTHOR

Dave Rolsky <autarch@urth.org>

=head1 COPYRIGHT AND LICENSE

This software is Copyright (c) 2015 by Dave Rolsky.

This is free software, licensed under:

  The Artistic License 2.0 (GPL Compatible)

=cut
