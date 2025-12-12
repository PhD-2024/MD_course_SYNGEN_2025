local pandoc = require 'pandoc'

local function safe_write(msg)
  pcall(function() io.stderr:write(msg) end)
end

local function clean_html_tags(s)
  if not s then return "" end
  s = s:gsub('<[^>]+>', '')
  s = s:gsub('&nbsp;', ' ')
  s = s:gsub('&lt;', '<'):gsub('&gt;', '>'):gsub('&amp;', '&')
  s = s:match('^%s*(.-)%s*$') or s
  return s
end

local function str_to_inlines(s)
  local inls = {}
  for word, sep in s:gmatch('(%S+)(%s*)') do
    table.insert(inls, pandoc.Str(word))
    if #sep > 0 then table.insert(inls, pandoc.Space()) end
  end
  return inls
end

function Pandoc(doc)
  safe_write('[details-debug] filter start\n')

  local path = '/home_w/dullinger/md_course_december_2025/MD_course/README.md'
  local f = io.open(path, 'r')
  local src = nil
  if f then
    src = f:read('*a')
    f:close()
    safe_write('[details-debug] scanned README.md on disk for headings\n')
  else
    safe_write('[details-debug] cannot open README.md on disk\n')
  end

  local headings = {}
  if src then
    for full, inner in src:gmatch('(<span[^>]-role="heading"[^>]->([%s%S]-)</span>)') do
      local cleaned = clean_html_tags(inner)
      if cleaned ~= '' then
        table.insert(headings, cleaned)
        safe_write(string.format('[details-debug] src heading found: %s\n', cleaned))
      end
    end
    safe_write(string.format('[details-debug] total src headings: %d\n', #headings))
  end

  if #headings == 0 then
    safe_write('[details-debug] no headings found in source; nothing to inject\n')
    return doc
  end

  local out_blocks = {}
  local h_idx = 1
  local i = 1
  local utils = pandoc.utils

  local function block_text_content(nb)
    if not nb then return nil end
    if nb.t == 'Para' or nb.t == 'Plain' then
      return utils.stringify(nb)
    elseif nb.t == 'RawBlock' and nb.format == 'html' then
      return clean_html_tags(nb.text)
    end
    return nil
  end

  while i <= #doc.blocks do
    local b = doc.blocks[i]
    -- If rendering to LaTeX/PDF, insert a page break before top-level sections (Header level 2)
    if FORMAT and tostring(FORMAT):match('latex') and b and b.t == 'Header' and b.level == 2 then
      table.insert(out_blocks, pandoc.RawBlock('latex', '\\clearpage'))
      safe_write('[details-debug] inserted \\clearpage before level-2 header\n')
    end
    if b.t == 'RawBlock' and b.format == 'html' and b.text:match('^%s*<details') then
      if h_idx <= #headings then
        local title = headings[h_idx]
        local inls = str_to_inlines(title)
        table.insert(out_blocks, pandoc.Header(3, inls))
        safe_write(string.format('[details-debug] injected Header #%d: %s\n', h_idx, title))
        h_idx = h_idx + 1
      else
        safe_write('[details-debug] no remaining src heading to inject for this <details>\n')
      end

      -- keep the original details block
      table.insert(out_blocks, b)

      -- inspect next block; if it contains the same heading text (or starts with it), skip it to avoid duplication
      local nb = doc.blocks[i+1]
      local next_text = block_text_content(nb)
      if next_text and #next_text > 0 then
        local title_lower = (headings[h_idx-1] or ''):lower()
        local next_lower = next_text:lower()
        if title_lower ~= '' and (next_lower == title_lower or next_lower:find(title_lower, 1, true)) then
          safe_write(string.format('[details-debug] skipped duplicate following block for heading: %s\n', title_lower))
          i = i + 2 -- skip the duplicate block
        else
          i = i + 1
        end
      else
        i = i + 1
      end
    else
      table.insert(out_blocks, b)
      i = i + 1
    end
  end

  safe_write(string.format('[details-debug] injected %d headers (requested %d)\n', h_idx-1, #headings))
  safe_write('[details-debug] filter end\n')

  return pandoc.Pandoc(out_blocks, doc.meta)
end