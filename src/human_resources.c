#include <ctype.h>
#include <dirent.h>
#include <errno.h>
#include <limits.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <zlib.h>
#include "human_resources.h"
#include "kommon.h"

#ifndef PATH_MAX
#define PATH_MAX 4096
#endif

typedef struct {
	const char *id;
	const char *rel_path;
	const char *basename;
} mb_hr_def_t;

typedef struct {
	uint64_t records, bases, checksum;
	char *line;
	size_t line_len, line_cap;
} mb_hr_stat_t;

typedef struct {
	const mb_hr_def_t *def;
	mb_hr_stat_t st;
	int found, build;
	char path[PATH_MAX];
} mb_hr_res_t;

typedef struct {
	mb_hr_res_t *res;
	int n_res;
} mb_hr_find_t;

static const mb_hr_def_t mb_hr_defs[] = {
	{ "mappability_150", "dna_pipeline/variants/mappability_150.38.bed.gz", "mappability_150.38.bed.gz" },
	{ "unmap_regions", "dna_pipeline/common/unmap_regions.38.tsv", "unmap_regions.38.tsv" },
	{ "sv_prep_blacklist", "dna_pipeline/sv/sv_prep_blacklist.38.bed", "sv_prep_blacklist.38.bed" },
	{ "hla", "other/lilac/hla.38.bed", "hla.38.bed" },
	{ "known_blacklist_germline", "dna_pipeline/variants/KnownBlacklist.germline.38.bed", "KnownBlacklist.germline.38.bed" }
};

static uint64_t mb_hr_fnv1a(uint64_t h, const unsigned char *buf, size_t len)
{
	size_t i;
	for (i = 0; i < len; ++i)
		h = (h ^ buf[i]) * 1099511628211ULL;
	return h;
}

static int mb_hr_name_build(const char *s)
{
	if (strstr(s, ".38.") || strstr(s, ".38/") || strstr(s, "_38.") || strstr(s, "GRCh38") || strstr(s, "grch38") || strstr(s, "hg38"))
		return 38;
	if (strstr(s, ".37.") || strstr(s, ".37/") || strstr(s, "_37.") || strstr(s, "GRCh37") || strstr(s, "grch37") || strstr(s, "hg19"))
		return 37;
	return 0;
}

static const char *mb_hr_build_name(int build)
{
	return build == 38? "GRCh38" : build == 37? "GRCh37" : "unknown";
}

static int mb_hr_is_chr(const char *name, const char *chr)
{
	if (strcmp(name, chr) == 0) return 1;
	if (strncmp(chr, "chr", 3) == 0 && strcmp(name, chr + 3) == 0) return 1;
	if (strncmp(name, "chr", 3) == 0 && strcmp(name + 3, chr) == 0) return 1;
	return 0;
}

static int mb_hr_index_build(const l2b_t *l2b)
{
	uint64_t i;
	int vote37 = 0, vote38 = 0;
	for (i = 0; i < l2b->n_ctg; ++i) {
		const char *n = l2b->ctg[i].name;
		uint64_t len = l2b->ctg[i].len;
		if (mb_hr_is_chr(n, "chr1")) {
			if (len == 248956422ULL) ++vote38;
			else if (len == 249250621ULL) ++vote37;
		} else if (mb_hr_is_chr(n, "chr2")) {
			if (len == 242193529ULL) ++vote38;
			else if (len == 243199373ULL) ++vote37;
		} else if (mb_hr_is_chr(n, "chrX")) {
			if (len == 156040895ULL) ++vote38;
			else if (len == 155270560ULL) ++vote37;
		} else if (mb_hr_is_chr(n, "chrY")) {
			if (len == 57227415ULL) ++vote38;
			else if (len == 59373566ULL) ++vote37;
		}
	}
	if (vote38 > 0 && vote37 == 0) return 38;
	if (vote37 > 0 && vote38 == 0) return 37;
	return 0;
}

static int mb_hr_parse_interval(const char *line, uint64_t *len)
{
	const char *p = line, *q;
	char *endp;
	uint64_t st, en;
	if (line[0] == 0 || line[0] == '#') return 0;
	while (*p && !isspace((unsigned char)*p)) ++p;
	while (*p && isspace((unsigned char)*p)) ++p;
	if (*p == 0) return 0;
	st = strtoull(p, &endp, 10);
	if (endp == p) return 0;
	q = endp;
	while (*q && isspace((unsigned char)*q)) ++q;
	if (*q == 0) return 0;
	en = strtoull(q, &endp, 10);
	if (endp == q || en < st) return 0;
	*len = en - st;
	return 1;
}

static void mb_hr_finish_line(mb_hr_stat_t *st)
{
	uint64_t len = 0;
	if (st->line_len > 0 && st->line[st->line_len - 1] == '\r')
		st->line[--st->line_len] = 0;
	if (st->line_len == 0) return;
	st->line[st->line_len] = 0;
	if (st->line[0] == '#') {
		st->line_len = 0;
		return;
	}
	++st->records;
	if (mb_hr_parse_interval(st->line, &len))
		st->bases += len;
	st->line_len = 0;
}

static void mb_hr_stat_update(mb_hr_stat_t *st, const unsigned char *buf, size_t len)
{
	size_t i;
	st->checksum = mb_hr_fnv1a(st->checksum, buf, len);
	for (i = 0; i < len; ++i) {
		unsigned char c = buf[i];
		if (c == '\n') {
			mb_hr_finish_line(st);
			continue;
		}
		if (st->line_len + 1 >= st->line_cap) {
			st->line_cap = st->line_cap? st->line_cap << 1 : 256;
			st->line = kom_realloc(char, st->line, st->line_cap);
		}
		st->line[st->line_len++] = c;
	}
}

static void mb_hr_stat_finish(mb_hr_stat_t *st)
{
	if (st->line_len > 0)
		mb_hr_finish_line(st);
}

static void mb_hr_stat_destroy(mb_hr_stat_t *st)
{
	free(st->line);
}

static int mb_hr_endswith(const char *s, const char *suffix)
{
	size_t l = strlen(s), m = strlen(suffix);
	return l >= m && strcmp(s + l - m, suffix) == 0;
}

static const char *mb_hr_basename(const char *path)
{
	const char *p = strrchr(path, '/');
	return p? p + 1 : path;
}

static int mb_hr_matches_def(const mb_hr_def_t *def, const char *path)
{
	const char *base = mb_hr_basename(path);
	char alt[PATH_MAX];
	const char *p;
	if (strcmp(base, def->basename) == 0) return 1;
	if (mb_hr_endswith(path, def->rel_path)) return 1;
	p = strstr(def->basename, ".38.");
	if (p) {
		size_t left = p - def->basename;
		snprintf(alt, sizeof(alt), "%.*s.37.%s", (int)left, def->basename, p + 4);
		if (strcmp(base, alt) == 0) return 1;
	}
	p = strstr(def->rel_path, ".38.");
	if (p) {
		size_t left = p - def->rel_path;
		snprintf(alt, sizeof(alt), "%.*s.37.%s", (int)left, def->rel_path, p + 4);
		if (mb_hr_endswith(path, alt)) return 1;
	}
	return 0;
}

static int mb_hr_resource_index(mb_hr_res_t *res, int n_res, const char *path)
{
	int i;
	for (i = 0; i < n_res; ++i)
		if (mb_hr_matches_def(res[i].def, path))
			return i;
	return -1;
}

static int mb_hr_read_plain_gz(const char *fn, mb_hr_stat_t *st)
{
	unsigned char buf[32768];
	gzFile fp;
	int n;
	st->checksum = 1469598103934665603ULL;
	fp = gzopen(fn, "rb");
	if (fp == 0) return -1;
	while ((n = gzread(fp, buf, sizeof(buf))) > 0)
		mb_hr_stat_update(st, buf, (size_t)n);
	if (gzclose(fp) != Z_OK) return -1;
	mb_hr_stat_finish(st);
	return 0;
}

static int mb_hr_read_nested_gz(const unsigned char *buf, size_t len, mb_hr_stat_t *st)
{
	z_stream zs;
	unsigned char out[32768];
	int ret;
	memset(&zs, 0, sizeof(zs));
	st->checksum = 1469598103934665603ULL;
	zs.next_in = (Bytef*)buf;
	zs.avail_in = (uInt)len;
	if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) return -1;
	do {
		zs.next_out = out;
		zs.avail_out = sizeof(out);
		ret = inflate(&zs, Z_NO_FLUSH);
		if (ret != Z_OK && ret != Z_STREAM_END) {
			inflateEnd(&zs);
			return -1;
		}
		mb_hr_stat_update(st, out, sizeof(out) - zs.avail_out);
	} while (ret != Z_STREAM_END);
	inflateEnd(&zs);
	mb_hr_stat_finish(st);
	return 0;
}

static int mb_hr_process_tar_member(gzFile tar, uint64_t size, int gzipped, mb_hr_stat_t *st)
{
	unsigned char buf[32768];
	uint64_t left = size;
	st->checksum = 1469598103934665603ULL;
	if (gzipped) {
		unsigned char *mem = 0;
		uint64_t off = 0, cap = size? size : 1;
		mem = kom_malloc(unsigned char, cap);
		while (left > 0) {
			unsigned int want = left > sizeof(buf)? sizeof(buf) : (unsigned int)left;
			int n = gzread(tar, buf, want);
			if (n <= 0) { free(mem); return -1; }
			if (off + (uint64_t)n > cap) {
				cap = off + n;
				mem = kom_realloc(unsigned char, mem, cap);
			}
			memcpy(mem + off, buf, n);
			off += n, left -= n;
		}
		if (mb_hr_read_nested_gz(mem, (size_t)off, st) < 0) {
			free(mem);
			return -1;
		}
		free(mem);
	} else {
		while (left > 0) {
			unsigned int want = left > sizeof(buf)? sizeof(buf) : (unsigned int)left;
			int n = gzread(tar, buf, want);
			if (n <= 0) return -1;
			mb_hr_stat_update(st, buf, (size_t)n);
			left -= n;
		}
		mb_hr_stat_finish(st);
	}
	return 0;
}

static int mb_hr_skip_gz(gzFile fp, uint64_t size)
{
	unsigned char buf[32768];
	while (size > 0) {
		unsigned int want = size > sizeof(buf)? sizeof(buf) : (unsigned int)size;
		int n = gzread(fp, buf, want);
		if (n <= 0) return -1;
		size -= n;
	}
	return 0;
}

static uint64_t mb_hr_tar_size(const unsigned char *hdr)
{
	uint64_t size = 0;
	int i;
	for (i = 124; i < 136 && hdr[i]; ++i)
		if (hdr[i] >= '0' && hdr[i] <= '7')
			size = (size << 3) + (hdr[i] - '0');
	return size;
}

static int mb_hr_tar_name(const unsigned char *hdr, char *name, size_t name_cap)
{
	char prefix[156], base[101];
	size_t lp, lb;
	memcpy(base, hdr, 100);
	base[100] = 0;
	memcpy(prefix, hdr + 345, 155);
	prefix[155] = 0;
	lp = strlen(prefix), lb = strlen(base);
	if (lb == 0) return -1;
	if (lp > 0) snprintf(name, name_cap, "%s/%s", prefix, base);
	else snprintf(name, name_cap, "%s", base);
	return 0;
}

static int mb_hr_import_tar(const char *fn, mb_hr_res_t *res, int n_res)
{
	unsigned char hdr[512];
	gzFile fp;
	int found = 0;
	fp = gzopen(fn, "rb");
	if (fp == 0) return -1;
	for (;;) {
		char name[PATH_MAX];
		uint64_t size, pad;
		int i, all_zero = 1, rid;
		int n = gzread(fp, hdr, 512);
		if (n == 0) break;
		if (n != 512) {
			fprintf(stderr, "ERROR: truncated tar header in '%s'\n", fn);
			gzclose(fp);
			return -1;
		}
		for (i = 0; i < 512; ++i)
			if (hdr[i] != 0) { all_zero = 0; break; }
		if (all_zero) break;
		if (mb_hr_tar_name(hdr, name, sizeof(name)) < 0) {
			fprintf(stderr, "ERROR: invalid tar member name in '%s'\n", fn);
			gzclose(fp);
			return -1;
		}
		size = mb_hr_tar_size(hdr);
		pad = (512 - (size & 511)) & 511;
		rid = (hdr[156] == 0 || hdr[156] == '0')? mb_hr_resource_index(res, n_res, name) : -1;
		if (rid >= 0 && !res[rid].found) {
			int gzipped = mb_hr_endswith(name, ".gz");
			if (mb_hr_process_tar_member(fp, size, gzipped, &res[rid].st) < 0) {
				fprintf(stderr, "ERROR: failed to read tar member '%s'\n", name);
				gzclose(fp);
				return -1;
			}
			res[rid].found = 1;
			res[rid].build = mb_hr_name_build(name);
			snprintf(res[rid].path, sizeof(res[rid].path), "%s", name);
			++found;
		} else if (mb_hr_skip_gz(fp, size) < 0) {
			fprintf(stderr, "ERROR: failed to skip tar member '%s'\n", name);
			gzclose(fp);
			return -1;
		}
		if (pad && mb_hr_skip_gz(fp, pad) < 0) {
			fprintf(stderr, "ERROR: failed to skip tar padding after '%s'\n", name);
			gzclose(fp);
			return -1;
		}
	}
	gzclose(fp);
	return found;
}

static int mb_hr_find_dir_cb(const char *path, const char *name, void *data)
{
	mb_hr_find_t *find = (mb_hr_find_t*)data;
	int rid = mb_hr_resource_index(find->res, find->n_res, name);
	if (rid < 0 || find->res[rid].found) return 0;
	if (mb_hr_read_plain_gz(path, &find->res[rid].st) < 0) return -1;
	find->res[rid].found = 1;
	find->res[rid].build = mb_hr_name_build(name);
	snprintf(find->res[rid].path, sizeof(find->res[rid].path), "%s", path);
	return 0;
}

static int mb_hr_walk_dir(const char *dir, int (*cb)(const char*, const char*, void*), void *data)
{
	DIR *dp;
	struct dirent *de;
	dp = opendir(dir);
	if (dp == 0) return -1;
	while ((de = readdir(dp)) != 0) {
		char path[PATH_MAX];
		struct stat st;
		if (strcmp(de->d_name, ".") == 0 || strcmp(de->d_name, "..") == 0) continue;
		snprintf(path, sizeof(path), "%s/%s", dir, de->d_name);
		if (stat(path, &st) != 0) { closedir(dp); return -1; }
		if (S_ISDIR(st.st_mode)) {
			if (mb_hr_walk_dir(path, cb, data) < 0) { closedir(dp); return -1; }
		} else if (S_ISREG(st.st_mode)) {
			if (cb(path, de->d_name, data) < 0) { closedir(dp); return -1; }
		}
	}
	closedir(dp);
	return 0;
}

static int mb_hr_import_dir(const char *dir, mb_hr_res_t *res, int n_res)
{
	int i, found = 0;
	mb_hr_find_t find;
	find.res = res, find.n_res = n_res;
	if (mb_hr_walk_dir(dir, mb_hr_find_dir_cb, &find) < 0) return -1;
	for (i = 0; i < n_res; ++i)
		if (res[i].found) ++found;
	return found;
}

static int mb_hr_resource_build(mb_hr_res_t *res, int n_res)
{
	int i, build = 0;
	for (i = 0; i < n_res; ++i) {
		if (!res[i].found || res[i].build == 0) continue;
		if (build == 0) build = res[i].build;
		else if (build != res[i].build) return -1;
	}
	return build;
}

static int mb_hr_write_manifest(const char *prefix, const char *src, mb_hr_res_t *res, int n_res, int index_build, int resource_build)
{
	char *fn;
	FILE *fp;
	int i;
	fn = kom_calloc(char, strlen(prefix) + 16);
	strcat(strcpy(fn, prefix), ".hmf.tsv");
	fp = fopen(fn, "wb");
	if (fp == 0) {
		free(fn);
		return -1;
	}
	fprintf(fp, "#minibwa_hmf_resources\t1\n");
	fprintf(fp, "#source\t%s\n", src);
	fprintf(fp, "#index_build\t%s\n", mb_hr_build_name(index_build));
	fprintf(fp, "#resource_build\t%s\n", mb_hr_build_name(resource_build));
	fprintf(fp, "resource\tbuild\trecords\tbases\tchecksum64\tpath\n");
	for (i = 0; i < n_res; ++i)
		fprintf(fp, "%s\t%s\t%lu\t%lu\t%016lx\t%s\n", res[i].def->id, mb_hr_build_name(res[i].build),
			(unsigned long)res[i].st.records, (unsigned long)res[i].st.bases,
			(unsigned long)res[i].st.checksum, res[i].path);
	fclose(fp);
	if (kom_verbose >= 3) fprintf(stderr, "[M::%s] wrote %s\n", __func__, fn);
	free(fn);
	return 0;
}

int mb_human_resources_import(const char *src, const char *prefix, const l2b_t *l2b)
{
	mb_hr_res_t res[sizeof(mb_hr_defs) / sizeof(mb_hr_defs[0])];
	struct stat st;
	int i, n_res = (int)(sizeof(mb_hr_defs) / sizeof(mb_hr_defs[0]));
	int found, index_build, resource_build;
	int ret = -1;
	memset(res, 0, sizeof(res));
	for (i = 0; i < n_res; ++i)
		res[i].def = &mb_hr_defs[i];
	if (stat(src, &st) != 0) {
		fprintf(stderr, "ERROR: failed to stat human resources '%s': %s\n", src, strerror(errno));
		return -1;
	}
	found = S_ISDIR(st.st_mode)? mb_hr_import_dir(src, res, n_res) : mb_hr_import_tar(src, res, n_res);
	if (found < 0) {
		fprintf(stderr, "ERROR: failed to read human resources '%s'\n", src);
		goto end_import;
	}
	if (found != n_res) {
		fprintf(stderr, "ERROR: found %d/%d HMF resources in '%s'\n", found, n_res, src);
		for (i = 0; i < n_res; ++i)
			if (!res[i].found) fprintf(stderr, "ERROR: missing %s\n", res[i].def->basename);
		goto end_import;
	}
	resource_build = mb_hr_resource_build(res, n_res);
	if (resource_build < 0) {
		fprintf(stderr, "ERROR: HMF resources mix genome builds\n");
		goto end_import;
	}
	index_build = mb_hr_index_build(l2b);
	if (index_build && resource_build && index_build != resource_build) {
		fprintf(stderr, "ERROR: HMF resource build %s does not match index build %s\n",
			mb_hr_build_name(resource_build), mb_hr_build_name(index_build));
		goto end_import;
	}
	if (kom_verbose >= 3) {
		for (i = 0; i < n_res; ++i)
			fprintf(stderr, "[M::%s] %s build=%s records=%lu bases=%lu checksum64=%016lx\n",
				__func__, res[i].def->id, mb_hr_build_name(res[i].build),
				(unsigned long)res[i].st.records, (unsigned long)res[i].st.bases,
				(unsigned long)res[i].st.checksum);
		if (index_build == 0 && resource_build != 0)
			fprintf(stderr, "[WARNING] HMF resources are %s but index genome build could not be inferred\n", mb_hr_build_name(resource_build));
	}
	if (mb_hr_write_manifest(prefix, src, res, n_res, index_build, resource_build) < 0) {
		fprintf(stderr, "ERROR: failed to write HMF resource sidecar for prefix '%s'\n", prefix);
		goto end_import;
	}
	ret = 0;
end_import:
	for (i = 0; i < n_res; ++i)
		mb_hr_stat_destroy(&res[i].st);
	return ret;
}

void mb_human_resources_warn_mismatch(const char *prefix, const l2b_t *l2b)
{
	char *fn, buf[256];
	FILE *fp;
	int index_build, resource_build = 0;
	fn = kom_calloc(char, strlen(prefix) + 16);
	strcat(strcpy(fn, prefix), ".hmf.tsv");
	fp = fopen(fn, "rb");
	if (fp == 0) {
		free(fn);
		return;
	}
	while (fgets(buf, sizeof(buf), fp) != 0) {
		if (strncmp(buf, "#resource_build\t", 16) == 0) {
			if (strstr(buf + 16, "GRCh38")) resource_build = 38;
			else if (strstr(buf + 16, "GRCh37")) resource_build = 37;
			break;
		}
	}
	fclose(fp);
	index_build = mb_hr_index_build(l2b);
	if (index_build && resource_build && index_build != resource_build && kom_verbose >= 2)
		fprintf(stderr, "[WARNING] HMF resource sidecar %s declares %s but index looks like %s\n",
			fn, mb_hr_build_name(resource_build), mb_hr_build_name(index_build));
	free(fn);
}
