# 🚂 Deploy to Railway Using CLI - Step by Step

## ⚠️ Important: Run from Your Local Machine

This environment is sandboxed and cannot access Railway's servers. You need to run these commands on **your local machine** where you have:
- ✅ Network access
- ✅ The repository cloned
- ✅ Node.js installed

---

## 🎯 Quick Deploy (Automated)

### One-Command Deployment

```bash
# From your local machine, in the repository directory:
./deploy-railway-cli.sh
```

This script handles everything automatically! Just answer a few prompts.

---

## 📋 Manual Step-by-Step (If You Prefer)

### Step 1: Prerequisites

**On Your Local Machine:**

```bash
# 1. Verify Node.js is installed
node -v  # Should be v18 or higher

# 2. Clone the repository (if not already done)
git clone https://github.com/MichaelCrowe11/crios-dr-crowe-coder.git
cd crios-dr-crowe-coder

# 3. Checkout the deployment branch
git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH
git pull origin claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH
```

---

### Step 2: Install Railway CLI

```bash
# Install globally via npm
npm i -g @railway/cli

# Verify installation
railway --version
# Should output: railway version x.x.x
```

---

### Step 3: Login to Railway

```bash
# This will open your browser for authentication
railway login

# Verify you're logged in
railway whoami
# Should show your Railway username/email
```

---

### Step 4: Initialize Railway Project

```bash
# Option A: Create new project
railway init
# Enter a project name (e.g., "crios-dr-crowe-coder")

# Option B: Link existing project
railway link
# Select from your existing Railway projects
```

---

### Step 5: Add Databases

```bash
# Add PostgreSQL (takes ~1-2 minutes to provision)
railway add --plugin postgresql

# Add Redis (takes ~1 minute to provision)
railway add --plugin redis

# Verify databases are provisioned
railway variables
# You should see DATABASE_URL, REDIS_URL, and other DB variables
```

---

### Step 6: Set Environment Variables

```bash
# Set your Anthropic API key (REQUIRED)
railway variables set ANTHROPIC_API_KEY=sk-ant-your_key_here

# Generate and set JWT secret
railway variables set JWT_SECRET=$(openssl rand -base64 32)

# Set production environment
railway variables set ENVIRONMENT=production
railway variables set DEBUG=false
railway variables set LOG_LEVEL=info
railway variables set PYTHONUNBUFFERED=1

# Verify all variables are set
railway variables
```

---

### Step 7: Deploy to Railway

```bash
# Deploy all services (this takes 10-15 minutes first time)
railway up

# Railway will:
# ✅ Detect all services from railway.toml
# ✅ Build Docker images for each service
# ✅ Deploy Backend, Frontend, Worker, Core API
# ✅ Connect to PostgreSQL and Redis
# ✅ Start health checks
```

**What's Happening:**
- Building Backend (5-7 min): Installing RDKit, Python deps
- Building Frontend (2-3 min): Building Next.js app
- Building Worker (3-5 min): Installing PyTorch, ML libraries
- **Total: ~10-15 minutes**

---

### Step 8: Monitor Deployment

```bash
# Check deployment status
railway status

# Watch backend logs
railway logs --service crios-backend --follow

# Watch frontend logs
railway logs --service crios-frontend --follow

# Watch worker logs
railway logs --service crios-worker --follow
```

---

### Step 9: Verify Deployment

```bash
# Get service URLs
railway status

# Test backend health (replace with your actual URL)
curl https://crios-backend-production-xxxx.up.railway.app/health

# Should return: {"status":"healthy"}

# Open frontend in browser
railway open crios-frontend
# Or visit the URL shown in `railway status`
```

---

## 🎯 Expected Output

After successful deployment, `railway status` should show:

```
Project: crios-dr-crowe-coder

Services:
  ✓ crios-backend     Deployed    https://crios-backend-xxx.railway.app
  ✓ crios-frontend    Deployed    https://crios-frontend-xxx.railway.app
  ✓ crios-worker      Deployed
  ✓ postgres          Running
  ✓ redis             Running

Environment: production
```

---

## 🛠️ Troubleshooting

### Issue: "railway: command not found"

```bash
# Make sure Railway CLI is installed
npm i -g @railway/cli

# If still not working, try with npx
npx @railway/cli login
```

### Issue: "Not logged in"

```bash
# Login again
railway login

# This should open your browser
# If browser doesn't open, copy the URL from terminal
```

### Issue: "Build failed - RDKit not found"

This shouldn't happen with our Dockerfiles, but if it does:
- Check `Dockerfile.api` has `librdkit-dev` in apt-get install
- View build logs: `railway logs --service crios-backend`

### Issue: "Frontend can't connect to backend"

```bash
# Get your backend URL from railway status
BACKEND_URL=$(railway status --json | jq -r '.services[] | select(.name=="crios-backend") | .url')

# Set it for frontend
railway variables set NEXT_PUBLIC_API_URL=$BACKEND_URL --service crios-frontend

# Restart frontend
railway restart --service crios-frontend
```

### Issue: "Database connection failed"

```bash
# Verify DATABASE_URL is set
railway variables get DATABASE_URL

# Should start with: postgresql://...

# Test connection
railway run --service crios-backend python -c "import psycopg2; print('DB OK')"
```

---

## 📊 Deployment Timeline

| Phase | Duration | What's Happening |
|-------|----------|------------------|
| CLI Install | 1 min | Installing Railway CLI |
| Login | 30 sec | Browser authentication |
| Project Setup | 1 min | Creating Railway project |
| Database Provisioning | 2-3 min | PostgreSQL + Redis |
| Env Vars | 1 min | Setting configuration |
| Backend Build | 5-7 min | RDKit, Python, FastAPI |
| Frontend Build | 2-3 min | Next.js compilation |
| Worker Build | 3-5 min | PyTorch, ML libraries |
| **First Deploy Total** | **~15 min** | ☕ Coffee time |
| Subsequent Deploys | 3-5 min | Cached Docker layers |

---

## 🎓 Post-Deployment Tasks

### 1. Set Up Custom Domain (Optional)

```bash
# Add your domain
railway domain add yourdomain.com --service crios-frontend

# Railway will provide DNS settings
# Add CNAME record to your DNS provider
```

### 2. Enable Monitoring

```bash
# Add Sentry for error tracking (optional)
railway variables set SENTRY_DSN=your_sentry_dsn

# View metrics in Railway dashboard
railway open
```

### 3. Run Database Migrations

```bash
# SSH into backend
railway run --service crios-backend bash

# Run migrations
python manage.py migrate  # or your migration command

# Exit
exit
```

### 4. Load Sample Data (Optional)

```bash
railway run --service crios-backend python scripts/load_sample_data.py
```

---

## 🔄 Continuous Deployment

### Set Up Auto-Deploy from GitHub

1. **Go to Railway Dashboard**: `railway open`
2. **Select Service** (e.g., crios-backend)
3. **Settings → Source**
4. **Connect GitHub Repository**
5. **Select Branch**: `claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH`
6. **Enable Auto-Deploy**

Now every push to the branch triggers automatic deployment! 🎉

---

## 📚 Useful Railway CLI Commands

```bash
# View all services
railway status

# View environment variables
railway variables

# Set environment variable
railway variables set KEY=value

# Delete environment variable
railway variables delete KEY

# View logs (live)
railway logs --service crios-backend --follow

# SSH into container
railway run --service crios-backend bash

# Restart service
railway restart --service crios-backend

# Open Railway dashboard
railway open

# Deploy (after making changes)
railway up

# Connect to database
railway connect postgres

# Run one-off command
railway run --service crios-backend python manage.py shell

# Delete service
railway service delete crios-worker

# Unlink project
railway unlink
```

---

## 💰 Cost Management

```bash
# View usage and costs
railway open
# Go to Settings → Usage

# Scale down for development
railway service scale crios-backend --replicas 1

# Stop services when not in use
railway down  # (if available in your CLI version)
```

---

## 🆘 Getting Help

### Railway Support
- **Docs**: https://docs.railway.app
- **Discord**: https://discord.gg/railway
- **Status**: https://status.railway.app

### CriOS Issues
- **GitHub Issues**: https://github.com/MichaelCrowe11/crios-dr-crowe-coder/issues
- **Documentation**: See `RAILWAY_DEPLOYMENT.md`

---

## ✅ Deployment Checklist

- [ ] Node.js 18+ installed locally
- [ ] Repository cloned and on correct branch
- [ ] Railway CLI installed (`npm i -g @railway/cli`)
- [ ] Logged in to Railway (`railway login`)
- [ ] Project created or linked (`railway init`)
- [ ] PostgreSQL added (`railway add --plugin postgresql`)
- [ ] Redis added (`railway add --plugin redis`)
- [ ] ANTHROPIC_API_KEY set (`railway variables set`)
- [ ] JWT_SECRET set (`railway variables set`)
- [ ] Deployment initiated (`railway up`)
- [ ] Services healthy (`railway status`)
- [ ] Backend responding (`curl /health`)
- [ ] Frontend accessible (visit URL)

---

## 🎉 You're Ready!

Run the automated script:

```bash
cd /path/to/crios-dr-crowe-coder
git checkout claude/hello-world-011CUj5N8vfNVVQ2Vb5gxXHH
./deploy-railway-cli.sh
```

Or follow the manual steps above for full control.

**Your AI drug discovery platform will be live in ~15 minutes! 🚀💊🔬**
