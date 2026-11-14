nextflow.enable.types = true

params {
  samplesheet: Path
  settings: Path
  annotations: Path
  assay: String = 'insitu'
  baitmap: Path?
  keep_rejected: Boolean = false
}

process PARSE {
  input:
  record(library: String, raw: Path)
  settings: Path
  annotations: Path

  output:
  record(library: library, parsed: file("${library}.pairs.gz"))

  script:
  """
  python3 ${projectDir}/bin/hic_parse.py --settings ${settings} --annotations ${annotations} \\
    --input ${raw} --output ${library}.pairs.gz
  """
}

process DEDUP {
  input:
  record(library: String, parsed: Path)
  settings: Path

  stage:
  stageAs parsed, 'parsed.pairs.gz'

  output:
  record(library: library, dedup: file("${library}.pairs.gz"))

  script:
  """
  python3 ${projectDir}/bin/hic_dedup.py --settings ${settings} --input parsed.pairs.gz --output ${library}.pairs.gz
  """
}

process FILTER {
  input:
  record(library: String, dedup: Path)
  settings: Path
  baits: Path?
  keep: Boolean

  stage:
  stageAs dedup, 'dedup.pairs.gz'

  output:
  record(
    library: library,
    filtered: file(keep ? "${library}.classes/${library}.pairs.gz" : "${library}.pairs.gz"),
    classes: file("${library}.classes", optional: true)
  )

  script:
  def bait_opt = baits ? "--baitmap ${baits}" : ''
  def class_opt = keep ? "--classes ${library}.classes" : ''
  def target = keep ? "${library}.classes/${library}.pairs.gz" : "${library}.pairs.gz"
  """
  mkdir -p \$(dirname ${target})
  python3 ${projectDir}/bin/hic_filter.py --settings ${settings} --input dedup.pairs.gz --output ${target} \\
    ${bait_opt} ${class_opt}
  """
}

process QC {
  input:
  record(library: String, parsed: Path, dedup: Path, filtered: Path)
  settings: Path

  stage:
  stageAs parsed, 'parsed.pairs.gz'
  stageAs dedup, 'dedup.pairs.gz'
  stageAs filtered, 'filtered.pairs.gz'

  output:
  record(
    library: library,
    qc: file("${library}.qc.json"),
    verdict: eval("cat ${library}.verdict")
  )

  script:
  """
  python3 ${projectDir}/bin/hic_qc.py --settings ${settings} --library ${library} \\
    --parsed parsed.pairs.gz --dedup dedup.pairs.gz --filtered filtered.pairs.gz \\
    --output ${library}.qc.json --verdict ${library}.verdict
  """
}

workflow {
  main:
  if( params.assay != 'insitu' && params.assay != 'capture' )
    error "assay must be insitu or capture"
  if( params.assay == 'capture' && !params.baitmap )
    error "assay capture needs a baitmap"
  bait = params.assay == 'capture' ? params.baitmap : null

  sheet = params.samplesheet
  samples = channel.fromList(sheet.splitCsv(header: true))
    .map { row -> record(library: row.library, condition: row.condition, raw: sheet.resolveSibling(row.raw)) }

  parsed = PARSE(samples.map { s -> record(library: s.library, raw: s.raw) }, params.settings, params.annotations)
  deduped = DEDUP(parsed, params.settings)
  filtered = FILTER(deduped, params.settings, bait, params.keep_rejected)
  qc = QC(parsed.join(deduped, by: 'library').join(filtered, by: 'library'), params.settings)

  released = samples
    .join(parsed, by: 'library')
    .join(deduped, by: 'library')
    .join(filtered, by: 'library')
    .join(qc, by: 'library')
    .map { r ->
      record(
        library: r.library,
        condition: r.condition,
        passed: r.verdict == 'pass',
        parsed: r.parsed,
        dedup: r.dedup,
        filtered: r.filtered,
        classes: r.classes,
        design: bait
      )
    }

  reports = samples
    .join(qc, by: 'library')
    .map { r -> record(library: r.library, condition: r.condition, passed: r.verdict == 'pass', qc: r.qc) }

  publish:
  pairs = released
  qc = reports
}

output {
  pairs {
    path { r ->
      r.parsed >> (r.passed ? 'pairs/parsed/' : null)
      r.dedup >> (r.passed ? 'pairs/dedup/' : null)
      r.filtered >> (r.passed ? "pairs/filtered/${r.library}.pairs.gz" : null)
      r.classes >> (r.passed ? 'classes/' : null)
      r.design >> 'design/'
    }
    index {
      path 'pairs.csv'
      header true
    }
  }

  qc {
    path { r -> r.qc >> 'qc/' }
    index {
      path 'qc.csv'
      header true
    }
  }
}
