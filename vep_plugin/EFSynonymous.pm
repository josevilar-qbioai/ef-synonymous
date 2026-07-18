=head1 NAME

EFSynonymous — Ensembl VEP plugin: score termodinámico glass-box de variantes
sinónimas (EF-Synonymous, σ + G>A).

=head1 SYNOPSIS

 vep -i input.vcf --plugin EFSynonymous
 vep -i input.vcf --plugin EFSynonymous,python=/usr/bin/python3

Añade tres campos a la salida de VEP para las variantes con consecuencia
C<synonymous_variant>:

  EF_sigma   ΔΔG de apilamiento con signo (Turner, kcal/mol) en ±10 nt
  EF_score   P(patogénica) del modelo glass-box σ + G>A
  EF_acmg    evidencia ACMG ilustrativa (PP3 / BP4 / —)

=head1 DESCRIPTION

Puente hacia el paquete C<ef-synonymous> (mismo motor que el prototipo web,
verificado a 1e-6). Para cada variante sinónima extrae el CDS del transcrito y la
posición/alelos a nivel de CDS, y delega el cálculo en la CLI Python
(C<ef-syn score --json>).

NOTA de rendimiento (prototipo): invoca Python por variante. Para uso a escala,
la vía de producción es una C<cache> precomputada (todas las sinónimas por
transcrito) o un endpoint REST; este plugin es la implementación de referencia
de la integración.

Requiere: C<pip install ef-synonymous> y que C<ef-syn> esté en el PATH (o pasar
C<python=...> apuntando al intérprete donde está instalado).

=cut

package EFSynonymous;

use strict;
use warnings;
use JSON::PP;

use base qw(Bio::EnsEMBL::Variation::Utils::BaseVepPlugin);

sub new {
  my $class = shift;
  my $self = $class->SUPER::new(@_);
  my $params = $self->params_hash();
  $self->{python} = $params->{python} || 'python3';
  # comprobación blanda de disponibilidad
  my $ok = `$self->{python} -c "import ef_synonymous" 2>&1`;
  warn "[EFSynonymous] paquete ef_synonymous no importable con $self->{python}; ".
       "instala con 'pip install ef-synonymous'\n" if $?;
  return $self;
}

sub feature_types { return ['Transcript']; }

sub get_header_info {
  return {
    EF_sigma => 'EF-Synonymous: signed stacking ΔΔG (Turner, kcal/mol), ±10 nt',
    EF_score => 'EF-Synonymous: P(pathogenic), glass-box σ+G>A model (RUO)',
    EF_acmg  => 'EF-Synonymous: illustrative ACMG code (PP3/BP4/—), uncalibrated',
  };
}

sub run {
  my ($self, $tva) = @_;

  # solo variantes sinónimas
  my %cons = map { $_->SO_term => 1 } @{ $tva->get_all_OverlapConsequences };
  return {} unless $cons{synonymous_variant};

  my $tv = $tva->transcript_variation;
  my $tr = $tva->transcript;

  my $cds     = $tr->translateable_seq;          # CDS de la hebra codificante
  my $cds_pos = $tv->cds_start;                   # 1-based
  return {} unless $cds && $cds_pos;

  # alelos a nivel de CDS (VEP ya los da en la hebra del transcrito)
  my ($ref, $alt) = ($tva->base_variation_feature_overlap->reference_allele->variation_feature_seq,
                     $tva->variation_feature_seq);
  return {} unless defined $ref && defined $alt && length($ref) == 1 && length($alt) == 1;

  # temp file con el CDS (evita CDS largos en la línea de comandos)
  my $tmp = "/tmp/efsyn_$$_" . $cds_pos . ".cds";
  open(my $fh, '>', $tmp) or return {};
  print $fh $cds; close($fh);

  my $cmd = sprintf('%s -m ef_synonymous.cli score --cds-file %s --pos %d --ref %s --alt %s --json 2>/dev/null',
                    $self->{python}, $tmp, $cds_pos, $ref, $alt);
  my $out = `$cmd`;
  unlink $tmp;
  return {} unless $out;

  my $d = eval { decode_json($out) };
  return {} unless $d;

  return {
    EF_sigma => sprintf('%.3f', $d->{features}{sigma_signed_kcal_mol}),
    EF_score => sprintf('%.3f', $d->{score_sigma}),
    EF_acmg  => $d->{acmg}{code},
  };
}

1;
