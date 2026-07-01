using Logging
using Dates
using Printf

# Industry-standard structured logger for K-SpecPart: every record is emitted as
# a single timestamped, level-tagged, component-tagged line, e.g.
#
#   [2026-06-30 13:20:01.123] [INFO ] [K-SpecPart] Vertices : 91976
#
# Verbose, per-candidate details are emitted at DEBUG level and hidden by
# default; set KSPECPART_LOG_LEVEL=debug to see them.
struct SpecPartLogger <: AbstractLogger
    min_level::LogLevel
    io::IO
end

Logging.min_enabled_level(logger::SpecPartLogger) = logger.min_level
Logging.shouldlog(logger::SpecPartLogger, level, _module, group, id) = true
Logging.catch_exceptions(logger::SpecPartLogger) = false

function _level_tag(level::LogLevel)
    level >= Logging.Error && return "ERROR"
    level >= Logging.Warn  && return "WARN "
    level >= Logging.Info  && return "INFO "
    return "DEBUG"
end

function Logging.handle_message(logger::SpecPartLogger, level, message, _module,
                                group, id, filepath, line; kwargs...)
    ts = Dates.format(now(), "yyyy-mm-dd HH:MM:SS.sss")
    tag = _level_tag(level)
    print(logger.io, "[", ts, "] [", tag, "] [K-SpecPart] ", message, "\n")
    for (k, v) in kwargs
        print(logger.io, "[", ts, "] [", tag, "] [K-SpecPart]     ", k, " = ", v, "\n")
    end
    flush(logger.io)
end

# Install the K-SpecPart logger as the global logger. Honors KSPECPART_LOG_LEVEL
# (info|debug). Safe to call repeatedly.
function setup_logging(io::IO = stderr)
    level = lowercase(get(ENV, "KSPECPART_LOG_LEVEL", "info")) == "debug" ?
            Logging.Debug : Logging.Info
    global_logger(SpecPartLogger(level, io))
    return nothing
end

const LOG_RULE = repeat("=", 60)
const LOG_THIN = repeat("-", 60)

# Emit "label : value" with an aligned label column for readable key/value logs.
log_kv(label::AbstractString, value) = @info @sprintf("%-13s: %s", label, string(value))
