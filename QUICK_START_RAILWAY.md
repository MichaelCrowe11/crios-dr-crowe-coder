# Quick Start: Deploy CriOS to Railway in 5 Minutes

The fastest way to get CriOS Dr. Crowe Coder running on Railway.

## Prerequisites (2 minutes)

1. **Sign up for Railway**: https://railway.app
2. **Install Railway CLI**:
   ```bash
   npm i -g @railway/cli
   ```
3. **Get Anthropic API Key**: https://console.anthropic.com

## Deploy with One Command (3 minutes)

```bash
# 1. Login to Railway
railway login

# 2. Run the deployment script
./deploy-railway.sh
```

That's it! The script will:
- ✅ Create Railway project
- ✅ Provision PostgreSQL database
- ✅ Provision Redis cache
- ✅ Deploy all 5 services
- ✅ Set up environment variables
- ✅ Run health checks

## Or Deploy Manually (5-10 minutes)

### Step 1: Login
```bash
railway login
```

### Step 2: Create Project
```bash
railway init
```

### Step 3: Add Databases
```bash
railway add --plugin postgresql
railway add --plugin redis
```

### Step 4: Set Environment Variables
```bash
# Required
railway variables set ANTHROPIC_API_KEY=sk-ant-your_key_here
railway variables set JWT_SECRET=$(openssl rand -base64 32)

# Production settings
railway variables set \
  ENVIRONMENT=production \
  DEBUG=false \
  LOG_LEVEL=info
```

### Step 5: Deploy Services

#### Backend API
```bash
railway up --service crios-backend
```

#### Frontend
```bash
railway up --service crios-frontend
```

#### Celery Worker
```bash
railway up --service crios-worker
```

### Step 6: Verify Deployment
```bash
# Check status
railway status

# View logs
railway logs --service crios-backend --follow
```

## Access Your Application

After deployment completes:

1. **Get URLs**:
   ```bash
   railway status
   ```

2. **Open Frontend**:
   ```bash
   railway open crios-frontend
   ```

3. **Test Backend**:
   ```bash
   curl https://your-backend-url.railway.app/health
   ```

## What Gets Deployed?

| Service | Description | Resources |
|---------|-------------|-----------|
| **PostgreSQL** | Database for compound storage | Railway-managed |
| **Redis** | Cache & message broker | Railway-managed |
| **Backend API** | FastAPI + RDKit chemistry engine | 1 vCPU, 512MB |
| **Frontend** | Next.js 14 dashboard | 1 vCPU, 512MB |
| **Celery Worker** | ML pipeline & background jobs | 1 vCPU, 512MB |

## Cost Estimate

**Hobby Plan**: ~$10-15/month
- 5 services running 24/7
- PostgreSQL + Redis included
- 100GB bandwidth

**Pro Plan**: ~$25-35/month
- Better performance
- Custom domains
- Team collaboration
- Automatic backups

## Troubleshooting

### Build Failed?
```bash
# Check logs
railway logs --service crios-backend

# Common issues:
# 1. Missing dependencies → Check Dockerfile
# 2. RDKit not found → Ensure librdkit-dev is installed
# 3. Python version → Must be 3.11+
```

### Can't Access Frontend?
```bash
# Ensure NEXT_PUBLIC_API_URL is set
railway variables get NEXT_PUBLIC_API_URL

# Should be: https://crios-backend.railway.app
# Update if needed:
railway variables set NEXT_PUBLIC_API_URL=https://your-backend-url.railway.app
railway up --service crios-frontend
```

### Database Connection Error?
```bash
# Verify DATABASE_URL
railway variables get DATABASE_URL

# Test connection
railway run --service crios-backend python -c "import psycopg2; print('OK')"
```

## Next Steps

1. ✅ **Configure Custom Domain**:
   ```bash
   railway domain add yourdomain.com --service crios-frontend
   ```

2. ✅ **Set Up Monitoring**:
   - Add Sentry DSN to environment variables
   - Enable Railway metrics

3. ✅ **Scale Services**:
   ```bash
   railway service scale crios-backend --replicas 3
   ```

4. ✅ **Enable Auto-Deploy**:
   - Connect GitHub repository in Railway dashboard
   - Push to deploy automatically

## Useful Commands

```bash
# View all services
railway status

# SSH into service
railway run --service crios-backend bash

# View environment variables
railway variables

# Restart service
railway restart --service crios-backend

# Delete service
railway service delete crios-worker

# Export environment variables
railway variables > env-backup.txt
```

## Support

- **Full Guide**: See `RAILWAY_DEPLOYMENT.md`
- **Railway Docs**: https://docs.railway.app
- **Discord**: https://discord.gg/railway

---

**Deployment Time**: 5 minutes with script, 10 minutes manual

**Ready to deploy? Run `./deploy-railway.sh` now! 🚀**
